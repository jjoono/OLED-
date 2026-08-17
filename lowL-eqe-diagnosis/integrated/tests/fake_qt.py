import sys, types

class _Signal:
    def __init__(self, *a): pass
    def __get__(self, obj, owner): return self
    def connect(self, *a, **k): pass
    def emit(self, *a, **k): pass

class QThread(object):
    def __init__(self, parent=None): pass

qtcore = types.ModuleType("PySide6.QtCore")
qtcore.QThread = QThread
qtcore.Signal = lambda *a, **k: _Signal()
qtcore.Slot = lambda *a, **k: (lambda f: f)
qtcore.QMutex = type("QMutex", (), {"lock": lambda s: None, "unlock": lambda s: None})
qtcore.QRecursiveMutex = qtcore.QMutex
pyside = types.ModuleType("PySide6"); pyside.QtCore = qtcore
sys.modules["PySide6"] = pyside
sys.modules["PySide6.QtCore"] = qtcore

cf = types.ModuleType("core_functions")
cf.log_message = lambda m: print("  [log]", m)
def save_file(df, path, header_lines, save_header=False):
    with open(path, "w") as f:
        f.write("\n".join(header_lines))
        df.to_csv(f, index=False, header=False, sep="\t")
cf.save_file = save_file
sys.modules["core_functions"] = cf
