import fake_qt  # PySide6/core_functions 스텁 (import 순서 중요)
import sys, os
sys.path.insert(0, "../jvl/src")
import numpy as np

from autotube_measurement import AutotubeMeasurement

rng = np.random.default_rng(0)

class FakeSource:
    voltage = 0.0
    low_light = None
    def as_voltage_source(self, c): pass
    def set_voltage(self, v): self.voltage = float(v)
    def read_current(self):
        return 0.0 if self.voltage <= 2 else 1e-8 * 10 ** (self.voltage - 2)
    def activate_output(self): pass
    def deactivate_output(self): pass
    def set_low_light_mode(self, on): self.low_light = on; print("  [SMU] low_light =", on)
    def pd_signal(self):
        return 0.0 if self.voltage <= 2 else 40e-6 * 10 ** (self.voltage - 3.0)

class FakeDMM:
    n = 1
    def __init__(self, src): self.src = src
    def reset(self): print("  [DMM] reset")
    def configure_burst(self, s, multimeter_range=10, nplc=1): self.n = int(s); print("  [DMM] burst n=%d range=%s" % (s, multimeter_range))
    def read_burst(self): return self.src.pd_signal() + 8.2e-3 + rng.normal(0, 15e-6, self.n)
    def measure_voltage(self, r=0): return float(self.read_burst()[0])

class FakeUno:
    def init_serial_connection(self): pass
    def trigger_relay(self, r): print("  [relay]", r)
    def close(self): pass

class Obj:
    def __getattr__(self, k): return Obj()
    def __call__(self, *a, **k): return None

class FakeParent:
    progressBar = Obj(); aw_start_measurement_pushButton = Obj()
    def plot_autotube_measurement(self, *a): pass
    def unselect_all_pixels(self): print("  [parent] unselect_all_pixels")

import time as _t
_sleep = _t.sleep
_t.sleep = lambda s: _sleep(min(s, 0.001))  # 시뮬이므로 대기 단축

def run_mode(accurate):
    print("\n===== accurate_eqe_mode =", accurate, "=====")
    params = {
        "min_voltage": 2.4, "max_voltage": 3.4, "changeover_voltage": 3.05,
        "low_voltage_step": 0.3, "high_voltage_step": 0.3,
        "scan_compliance": 1050, "auto_spectrum": False,
        "photodiode_saturation": 10.0,
        "accurate_eqe_mode": accurate, "accurate_eqe_precision": 5.0,
    }
    setup = {"folder_path": "./out_", "batch_name": "smoke", "device_number": 1,
             "top_emitting": False}
    src = FakeSource()
    m = AutotubeMeasurement(src, FakeDMM(src), FakeUno(), params, setup,
                            0.001, [1, 2], 45.0, FakeParent())
    m.run()
    path = "./out_" + __import__("datetime").date.today().strftime("%Y-%m-%d_") + "smoke_d1_p2_jvl.csv"
    print("--- saved file (pixel 2) ---")
    print(open(path).read())
    assert src.low_light is (False if accurate else None)

run_mode(False)
run_mode(True)
print("\nSMOKE TEST PASSED")
