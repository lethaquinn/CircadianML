import sys
import os
sys.path.append('src')

def test_imports():
    try:
        from utils.data_loader import CircadianDataLoader
        from rhythm_detection.detector import CircadianRhythmDetector
        print('All imports successful')
        return True
    except ImportError as e:
        print(f'Import error: {e}')
        return False

if __name__ == '__main__':
    test_imports()
