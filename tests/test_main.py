from src import Varcall

# WARN: Valid only when samtools is focused on startup
def test_samtoolstab(snap_compare):
    assert snap_compare("../src/main.py")

def test_hometab(snap_compare):
    assert snap_compare("../src/main.py", press=["tab", "left"])

def test_bcftools(snap_compare):
    assert snap_compare("../src/main.py", press=["tab", "right"])

def test_helptab(snap_compare):
    assert snap_compare("../src/main.py", press=["tab", "right", "right"])

if __name__ == "__main__":
    test_samtoolstab(Varcall)
    test_hometab(Varcall)
