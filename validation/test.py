import numpy as np
from numpy.ma.core import arccos

def find_cosine(v1, v2):
    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

if __name__ == '__main__':
    # v1 = np.array([0.91276294686117521, 1.5809811135721659, 0])
    # v2 = np.array([0.91276294686117521, -1.5809811135721659, 0])
    # v3 = np.array([-1.8255675806274383, 0, 0])

    v1 = np.array([-0.0163449510307766, 0, 0])
    v2 = np.array([0.0081722888965992946, -0.014155082044696624, 0])
    v3 = np.array([0.0081722888965992946, 0.014155082044696624, 0])

    print(np.linalg.norm(v1))
    print(np.linalg.norm(v2))
    print(np.linalg.norm(v3))

    cos1 = find_cosine(v1, v3)
    cos2 = find_cosine(v2, v3)
    cos3 = find_cosine(v2, v1)

    print(cos1)
    print(cos2)
    print(cos3)

    print(arccos(cos1) * 180 / np.pi)
    print(arccos(cos2) * 180 / np.pi)
    print(arccos(cos3) * 180 / np.pi)

