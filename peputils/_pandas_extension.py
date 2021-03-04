import pandas as pd
import hashlib
import statsmodels.stats.multitest as multitest
import scipy as sp
import scipy.stats
import numpy as np
from collections import namedtuple
import warnings

from tqdm import tqdm

#  from .annotations import ProteinId, PeptideId, PeptideVariantId
#  from tqdm import tqdm

#  from modlamp.sequences import Helices
#  from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor


################################################################################
# namedtuples
################################################################################
TTestResult = namedtuple("TTestResult", ("pvalue", "statistic", "ratio"))
PeptideVariantId = namedtuple('PeptideVariantId',
                              ('protein_id', 'start', 'stop', 'mod_seq', 'origin'))
PeptideId = namedtuple('PeptideId', ('protein_id', 'start', 'stop'))
ProteinId = namedtuple('ProteinId', ('protein_id',))

 
################################################################################
# Pandas peptidomics Excensions
################################################################################
@pd.api.extensions.register_index_accessor("peptidomics")
class PandasPeptidomicsIndex:
    def __init__(self, index):
        self.index = index

        id_classes = (PeptideVariantId, PeptideId, ProteinId)
        id_levels = ("peptide_variants", "peptides", "proteins")
        for id_class, level in zip(id_classes, id_levels):
            if id_class._fields == self.index.names[1:]:
                self.id_class = id_class
                self.level = level
                break
        else:  # no break
            raise ValueError("index does not seem to be 'proteins', 'peptides' or "
                             "peptide_variants")

    def to_proteins(self):
        if self.level == 'peptides':
            drop = ("start", "stop")
        elif self.level == "peptide_variants":
            drop = ("start", "stop", "mod_seq", "origin")
        else:
            raise ValueError("index does not seem to be 'peptides' or 'peptide_variants'")
        return self._drop_index(drop)

    def to_peptides(self):
        if self.level != 'peptide_variants':
            raise ValueError("can only craete peptide index from a peptide_variant index")
        return self._drop_index(("mod_seq", "origin"))

    def _drop_index(self, drop):
        index = self.index.copy()
        for level in drop:
            index = index.droplevel(level)
        return index.unique()

    def iter_index(self):
        for row in self.index:
            yield row[0], self.id_class(*row[1:])

    def create_annotation_df(self):
        annotation_columns = pd.MultiIndex.from_tuples((), names=['annotation', 'group'])
        return pd.DataFrame(columns=annotation_columns, index=self.index)


@pd.api.extensions.register_series_accessor("peptidomics")
class PandasPeptidomicsSerie:
    def __init__(self, series):
        self.series = series

    def fdr(self, alpha=0.05, method='fdr_bh'):
        if method != "fdr_bh":
            raise NotImplementedError("TODO!!!")

        mask = self.series == self.series
        corrected = np.full(self.series.shape, np.nan)
        corrected[mask] = multitest.multipletests(self.series[mask], method=method, alpha=alpha)[1]
        return pd.Series(corrected, index=self.series.index, name=("FDR", self.series.name[1]))


@pd.api.extensions.register_dataframe_accessor("peptidomics")
class PandasPeptidomics:
    rename_header = {"end": "stop", "begin": "start", "prot_acc": "protein_id",
                     "pep_start": "start", "pep_end": "stop", "pep_mod_seq": "mod_seq",
                     "pep_seq": "seq"}

    def __init__(self, df):
        self.df = df

    def create_annotation(self):
        annotation_columns = pd.MultiIndex.from_tuples((), names=['annotation', 'group'])
        return pd.DataFrame(columns=annotation_columns, index=self.df.index)

    def add_annotation(self, annotation_series):
        self.df[annotation_series.name] = annotation_series

    @classmethod
    def load_ppv_file(cls, ppv_file, campaign_id, index=None):
        print('Load ppv file')
        ppv = pd.read_csv(ppv_file, sep="\t")
        ppv.rename(columns=cls.rename_header, inplace=True)
        ppv["campaign_id"] = campaign_id
        ppv["origin"] = "collapsed"
        if index is None:
            index = ['campaign_id'] + list(PeptideVariantId._fields)
        return pd.pivot_table(ppv, index=index, values="score")

    def _ppv_pval_combine(self, subset, compare, pvalue_name="P-Value", ratio_name="Ratio"):
        gt = subset[subset[ratio_name, compare] >= 0][pvalue_name, compare] / 2
        lt = subset[subset[ratio_name, compare] <= 0][pvalue_name, compare] / 2
        gt_flip = 1 - gt
        lt_flip = 1 - lt
        pval_gt = np.append(gt.values, lt_flip.values)
        pval_lt = np.append(lt.values, gt_flip.values)
        p1 = sp.stats.combine_pvalues(pval_gt)[1]
        p2 = sp.stats.combine_pvalues(pval_lt)[1]
        pval = min(p1, p2) * 2
        return pval

    def ppv_pval(self, df_observed, compare, method="fisher",
                 annotation_name="P-Value", progresbar=True):
        if method != "fisher":
            raise NotImplementedError("Only fishers method is implemented")
        _iter = (self._ppv_pval_combine(subset, compare)
                 for subset in self._iter_ppvs(df_observed))
        if progresbar:
            _iter = tqdm(_iter, "Collapsing p-values from UPF", total=self.df.shape[0])
        return pd.Series(list(_iter), index=self.df.index, name=(annotation_name, compare))

    def _ppv_intensity_combine(self, subset, group):
        return pd.Series(subset[group].values.flatten()).mean()

    def ppv_intensity(self, df_observed, group, annotation_name="Intensity", progresbar=True):
        _iter = (self._ppv_intensity_combine(subset, group)
                 for subset in self._iter_ppvs(df_observed))
        if progresbar:
            _iter = tqdm(_iter, "Collapsing intensities from UPF", total=self.df.shape[0])

        return pd.Series(list(_iter), index=self.df.index, name=(annotation_name, group))

    def _iter_ppvs(self, df):
        for index, ppv_data in self.df.iterrows():
            campaign_id, protein_id, start, stop, mod_seq, origin = index
            subset = df.query('protein_id == "{}" and start >= {} and stop <= {}'.format(
                              protein_id, start, stop))
            yield subset

    def ppv_ratio(self, df_observed, compare, annotation_name="Ratio", progresbar=True):
        _iter = (subset[annotation_name, compare].median()
                 for subset in self._iter_ppvs(df_observed))
        if progresbar:
            _iter = tqdm(_iter, "Annotating ratio from UPF", total=self.df.shape[0])
        return pd.Series(list(_iter), index=self.df.index, name=(annotation_name, compare))

    @classmethod
    def load_upf_meta(cls, upf_file, meta_file, campaign_id, upf_file_sep="\t", meta_file_sep="\t",
                      left_on='rs_acc', right_on='accno', pivot=True, add_meta_defaults=True,
                      **kwargs):
        # TODO: campaing ID should be infered from dataset.meta
        #  print('Load, merge and create peptide id')
        upf = pd.read_csv(upf_file, sep=upf_file_sep)
        upf.rename(columns=cls.rename_header, inplace=True)
        meta = cls._load_meta_from_file(meta_file, meta_file_sep, add_meta_defaults)
        upf_sample = pd.merge(meta, upf, left_on=left_on, right_on=right_on)
        upf_sample["campaign_id"] = campaign_id
        upf_sample["origin"] = "observed"
        if pivot:
            return cls.ms_pivot_table(upf_sample, **kwargs)
        return upf_sample

    @classmethod
    def _load_meta_from_file(cls, meta_file, sep, add_defaults=True):
        meta = pd.read_csv(meta_file, sep=sep)
        if not add_defaults:
            return meta
        if 'subject' not in meta.columns or meta['subject'].dropna().shape[0] == 0:
            # if there is no subject info assume that each sample is a unique subject
            meta["subject"] = meta["ps_acc"]
        if "qc" not in meta.columns or meta['qc'].dropna().shape[0] == 0:
            # if there are no qc info, all samples were probbably fine :D
            meta["qc"] = "OK"
        return meta

    @classmethod
    def ms_pivot_table(cls, df, values='intensity',
                       index=['campaign_id', 'protein_id', 'start', 'stop', 'mod_seq', 'origin'],
                       columns=['subject', 'group', 'qc']):
        return pd.pivot_table(df, values=values, index=index, columns=columns)

    def normalize(self, full_rank='auto', full_rank_cutoff=100):
        raw_data = self.df.replace(0, np.nan)
        log10_data = np.log10(raw_data)
        log10_data_full_rank = log10_data.dropna()
        rank = log10_data_full_rank.shape[0]
        if full_rank == 'auto':
            full_rank = rank >= full_rank_cutoff
        if full_rank:
            if rank < full_rank_cutoff:
                msg = ("The full rank of the normalization matrix is {}, consider "
                       "consider running df.peptidomics.normalize(full_rank=False) instead")
                warnings.warn(msg.format(rank))
            median_intensities = log10_data_full_rank.median()
        else:
            median_intensities = log10_data.median()
        norm_target = median_intensities.mean()
        normalization_factors = norm_target - median_intensities
        #  ax = full_rank_median_intensities.plot(kind="hist",
        #                                         title="Median Intensity across samples")
        return 10 ** (log10_data + normalization_factors)

    @classmethod
    def _ttest_return(cls, t_test, index, ratios):
        pvalue = pd.Series(t_test.pvalue.data, index=index, name="P-Value")
        statistic = pd.Series(t_test.statistic.data, index=index, name="T-Statistic")
        ratio = pd.Series(ratios, index=index, name="Ratio")
        return TTestResult(pvalue, statistic, ratio)

    def _paired_ttest(self, c1, c2, n_min):
        counts = (((c1 == c1).astype(int) + (c2 == c2).astype(int)) == 2).sum(axis=1)
        mask = n_min <= counts
        t_test = sp.stats.ttest_rel(c1[mask], c2[mask], axis=1, nan_policy='omit')
        ratios = (c1[mask] - c2[mask]).mean(axis=1)
        return self._ttest_return(t_test, self.df.index[mask], ratios)

    def _unpaired_ttest(self, c1, c2, n_min, equal_var):
        _data = {'c1': (c1 == c1).sum(axis=1), 'c2': (c2 == c2).sum(axis=1)}
        counts = pd.DataFrame(_data).min(axis=1)
        mask = n_min <= counts
        t_test = sp.stats.ttest_ind(c1[mask], c2[mask], axis=1, nan_policy='omit',
                                    equal_var=equal_var)
        ratios = c1[mask].mean(axis=1) - c2[mask].mean(axis=1)
        return self._ttest_return(t_test, self.df.index[mask], ratios)

    def calculate_t_test(self, cond1, cond2, ttype='unpaired', level='group', n_min=2,
                         equal_var=False):
        c1 = self.df.xs(cond1, level=level, axis=1)
        c2 = self.df.xs(cond2, level=level, axis=1)

        if ttype == 'paired':
            return self._paired_ttest(c1, c2, n_min)
        elif ttype == 'unpaired':
            return self._unpaired_ttest(c1, c2, n_min, equal_var)
        raise ValueError("ttype must be either 'paired' or 'unpaired' not {}".format(ttype))

    @classmethod
    def peptide_id_label(cls, acc, start, end, mod_seq):
        md5sum = hashlib.md5((mod_seq).encode('utf-8')).hexdigest().upper()
        return "%s:%s:%s:%s" % (acc, start, end, md5sum)

    def is_annotation(self):
        return self.df.columns.names == ['annotation', 'group']

    def iterrows(self):
        id_class = self.df.index.peptidomics.id_class
        for id_tuple, data in self.df.iterrows():
            yield id_class(*id_tuple[1:]), data

        #  for id_tuple, data in self.
