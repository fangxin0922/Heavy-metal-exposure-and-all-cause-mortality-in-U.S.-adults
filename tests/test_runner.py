"""Validate runner scheduling using mocked R processes; no numerical analyses."""
from contextlib import redirect_stdout
from importlib.util import module_from_spec, spec_from_file_location
from io import StringIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from threading import Event, Lock
from time import sleep
from types import SimpleNamespace
from unittest import TestCase, main
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parent.parent
SPEC = spec_from_file_location(
    'runner', REPOSITORY / 'metal_survival_analysis' / 'run_analysis.py')
runner = module_from_spec(SPEC)
sys.dont_write_bytecode = True
SPEC.loader.exec_module(runner)
EXPECTED_PRIMARY = ['primary_Pb', 'primary_Hg', 'egfr60_Cd', 'egfr60_As']
EXPECTED_ALL = [
    'primary_Pb', 'primary_Hg', 'egfr60_Cd', 'egfr60_As',
    'primary_Cd', 'primary_As',
    'no_bmi_Pb', 'no_bmi_Hg', 'no_bmi_egfr60_Cd', 'no_bmi_egfr60_As',
    'age55plus_Pb', 'age55plus_Hg', 'age55plus_egfr60_Cd', 'age55plus_egfr60_As',
    'portnoy_Pb', 'portnoy_Hg', 'portnoy_egfr60_Cd', 'portnoy_egfr60_As',
    'joint_egfr60_Pb', 'joint_egfr60_Cd', 'epa_dha_Hg', 'rcs_Pb', 'rcs_egfr60_Cd']


class RProcesses:
    def __init__(self, scope='all', failure=None):
        self.scope = scope
        self.failure = failure
        self.calls = []
        self.prepared = set()
        self.lock = Lock()
        self.fit_started = Event()
        self.secondary_started = Event()
        self.active_fits = 0
        self.peak_fits = 0
        self.overlap = False

    def __call__(self, command, **kwargs):
        assert command[:2] == ['Rscript', '--vanilla']
        script, arguments = Path(command[2]).name, command[3:]
        with self.lock:
            self.calls.append((script, arguments))
        kwargs['stdout'].write('Mock R process\n')
        if script == 'impute_covariates.R':
            scenario = arguments[0]
            if scenario not in EXPECTED_PRIMARY:
                with self.lock:
                    required = set(EXPECTED_PRIMARY) - {self.failure}
                    assert required <= self.prepared, 'Sensitivity preparation started early'
                    dependency = runner.PREPARATION_DEPENDENCIES.get(scenario)
                    assert dependency is None or dependency in self.prepared
                    self.overlap |= self.active_fits > 0
                self.secondary_started.set()
            elif scenario == 'egfr60_As' and self.scope == 'all':
                assert self.fit_started.wait(5), 'Primary fits were not scheduled promptly'
            if scenario == self.failure:
                return SimpleNamespace(returncode=1)
            with self.lock:
                self.prepared.add(scenario)
        elif script == 'fit_cqr.R':
            with self.lock:
                assert arguments[0] in self.prepared, 'Fit started before preparation finished'
                self.active_fits += 1
                self.peak_fits = max(self.peak_fits, self.active_fits)
            self.fit_started.set()
            if self.scope == 'all':
                assert self.secondary_started.wait(5), 'Fits blocked sensitivity preparation'
            sleep(0.003)
            with self.lock:
                self.active_fits -= 1
        return SimpleNamespace(returncode=0)


class RunnerTests(TestCase):
    def run_mocked(self, process, directory, *arguments):
        run = Path(directory) / 'metal_survival_analysis'
        run.mkdir(exist_ok=True)
        with patch.object(runner, 'ROOT', Path(directory)), \
                patch.object(runner, 'RUN', run), \
                patch.object(runner.subprocess, 'run', process), \
                patch('sys.argv', ['run_analysis.py', '--workers', '4',
                                   '--prepare-workers', '4', *arguments]), \
                redirect_stdout(StringIO()):
            runner.main()

    def test_all_scenarios_order_and_fit_concurrency(self):
        self.assertEqual(runner.ALL, EXPECTED_ALL)
        self.assertEqual(len(set(runner.ALL)), 23)
        process = RProcesses()
        with TemporaryDirectory() as directory:
            self.run_mocked(process, directory, '--points-only')
        preparations = [args[0] for script, args in process.calls
                        if script == 'impute_covariates.R']
        self.assertEqual(set(preparations[:4]), set(EXPECTED_PRIMARY))
        self.assertCountEqual(preparations, EXPECTED_ALL)
        fits = [args for script, args in process.calls if script == 'fit_cqr.R']
        self.assertEqual(len(fits), 23 * 5)
        for scenario in EXPECTED_ALL:
            self.assertCountEqual([args[1] for args in fits if args[0] == scenario],
                                  ['1', '2', '3', '4', '5'])
        self.assertTrue(all(args[2] == '0' for args in fits))
        self.assertGreater(process.peak_fits, 1)
        self.assertTrue(process.overlap)
        self.assertIn(('summarize_results.R', EXPECTED_ALL), process.calls)
        self.assertIn(('baseline_characteristics.R', []), process.calls)
        self.assertIn(('participant_flow.R', []), process.calls)
        self.assertIn(('supplementary_descriptives.R', []), process.calls)
        self.assertEqual([script for script, _ in process.calls[-4:]],
                         ['summarize_results.R', 'spline_predictions.R',
                          'survey_cox.R', 'draw_figures.R'])

    def test_primary_scope(self):
        process = RProcesses(scope='primary')
        with TemporaryDirectory() as directory:
            self.run_mocked(process, directory, '--scope', 'primary')
        self.assertEqual(runner.PRIMARY, EXPECTED_PRIMARY)
        self.assertEqual(process.prepared, set(EXPECTED_PRIMARY))
        fits = [args for script, args in process.calls if script == 'fit_cqr.R']
        self.assertEqual(len(fits), 20)
        self.assertTrue(all(args[2] == '400' for args in fits))
        self.assertIn(('summarize_results.R', EXPECTED_PRIMARY), process.calls)
        self.assertIn(('baseline_characteristics.R', []), process.calls)
        self.assertIn(('participant_flow.R', []), process.calls)
        self.assertIn(('supplementary_descriptives.R', []), process.calls)

    def test_failed_renal_primary_blocks_only_its_dependencies(self):
        for failed in ['egfr60_Cd', 'egfr60_As']:
            with self.subTest(failed=failed), TemporaryDirectory() as directory:
                process = RProcesses(failure=failed)
                with self.assertRaises(RuntimeError) as raised:
                    self.run_mocked(process, directory)
                blocked = {scenario for scenario, dependency
                           in runner.PREPARATION_DEPENDENCIES.items() if dependency == failed}
                attempted = {args[0] for script, args in process.calls
                             if script == 'impute_covariates.R'}
                self.assertEqual(attempted, set(EXPECTED_ALL) - blocked)
                fitted = {args[0] for script, args in process.calls if script == 'fit_cqr.R'}
                self.assertEqual(fitted, set(EXPECTED_ALL) - blocked - {failed})
                for scenario in blocked:
                    self.assertIn(f'{scenario} skipped: required preparation {failed} failed.',
                                  str(raised.exception))
                self.assertFalse(any(script == 'summarize_results.R'
                                     for script, _ in process.calls))

    def test_nonempty_output_safeguard(self):
        process = RProcesses(scope='primary')
        with TemporaryDirectory() as directory:
            self.run_mocked(process, directory, '--scope', 'primary')
            call_count = len(process.calls)
            with self.assertRaisesRegex(SystemExit, 'Use empty output folders'):
                self.run_mocked(process, directory, '--scope', 'primary')
            self.assertEqual(len(process.calls), call_count)


if __name__ == '__main__':
    main(verbosity=2)
