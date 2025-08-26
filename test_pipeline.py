import sys
sys.path.append('src')

from utils.data_loader import CircadianDataLoader
from rhythm_detection.detector import CircadianRhythmDetector

def test_pipeline():
    print('Testing CircadianML Pipeline...')
    loader = CircadianDataLoader()
    data = loader.generate_synthetic_data(n_genes=50)
    processed = loader.preprocess(data)
    
    detector = CircadianRhythmDetector()
    results = detector.detect_all(processed, loader.time_points)
    
    print(f'Detected {results["consensus_rhythmic"].sum()} rhythmic genes')
    return results

if __name__ == '__main__':
    test_pipeline()
    print('Pipeline test completed!')
