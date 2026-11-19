import importlib.util
from pathlib import Path
import unittest

path=Path(__file__).resolve().parents[1]/'validation/compare_results.py'
spec=importlib.util.spec_from_file_location('compare_results',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class ResultComparisonTests(unittest.TestCase):
    def row(self,**changes):
        row=dict(scenario='egfr60_Cd',imputation='1',term='exposure',tau='0.1',n='7061',
                 deaths='1510',estimate='-1.7',lower='-2.5',upper='-1.0')
        row.update(changes);return row

    def test_small_change_preserves_agreement(self):
        result=module.compare('main',[self.row()],[self.row(estimate='-1.71')])[0]
        self.assertFalse(result['magnitude_flag']);self.assertFalse(result['sign_changed'])

    def test_unavailable_is_not_a_matching_zero(self):
        result=module.compare('main',[self.row()],[self.row(estimate='',lower='',upper='')])[0]
        self.assertFalse(result['same_availability']);self.assertFalse(result['paired_finite'])
        self.assertTrue(result['ci_availability_changed'])

    def test_inference_change_is_flagged_even_with_same_point(self):
        result=module.compare('main',[self.row()],[self.row(upper='.01')])[0]
        self.assertTrue(result['ci_zero_classification_changed'])

    def test_missing_cells_and_changed_cohorts_fail(self):
        for rows in [[],[self.row(n='7060')],[self.row(),self.row()]]:
            with self.assertRaises(ValueError):module.compare('main',[self.row()],rows)

    def test_large_opposite_effect_is_flagged(self):
        result=module.compare('main',[self.row()],[self.row(estimate='1.7')])[0]
        self.assertTrue(result['magnitude_flag']);self.assertTrue(result['sign_changed'])

    def test_biomarker_interval_availability_and_replicate_count(self):
        reference=dict(pair='synthetic_marker',tau='.1',n='100',deaths='20',
          biomarker_low='1',biomarker_high='3',valid_replicates='399',p_value='',q_primary='',
          ic_years='.2',ic_lower_years='',ic_upper_years='',
          metal_low_years='-.5',metal_low_lower='',metal_low_upper='',
          metal_high_years='-.3',metal_high_lower='',metal_high_upper='')
        fresh=dict(reference,valid_replicates='400',ic_lower_years='-.1',ic_upper_years='.5',p_value='.2',q_primary='1')
        rows=module.compare('biomarker',[reference],[fresh])
        self.assertTrue(rows[0]['same_availability'])
        self.assertTrue(rows[0]['ci_availability_changed'])
        self.assertEqual(rows[0]['valid_replicates_difference'],1)
        self.assertTrue(rows[0]['p_availability_changed'])
        self.assertTrue(rows[0]['bh_availability_changed'])
        self.assertEqual(rows[0]['p_significance_changed'],'')
        self.assertEqual([r['valid_replicates_difference'] for r in rows[1:]],['',''])
        with self.assertRaisesRegex(ValueError,'reference level changed'):
            module.compare('biomarker',[reference],[dict(fresh,biomarker_high='4')])

if __name__=='__main__':unittest.main()
