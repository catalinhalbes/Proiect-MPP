import sys, os
import struct
import numpy as np

def loadmat(filename: str) -> np.array:
    with open(filename, 'rb') as f:
        dim_x, dim_y, dim_z = struct.unpack('<QQQ', f.read(24))
        data = np.frombuffer(f.read(), dtype=np.float64).reshape(dim_x, dim_y, dim_z)
        return data

def writemat(filename: str, mat: np.array) -> None:
    with open(filename, 'wb') as f:
        f.write(struct.pack("<QQ", mat.shape[0], mat.shape[1]))
        f.write(mat.tobytes())

def mat_MAPE(mat1_file, mat2_file) -> float:
    mat1 = loadmat(mat1_file)
    mat2 = loadmat(mat2_file)

    mask = mat1 != 0
    mat1 = mat1[mask]
    mat2 = mat2[mask]
    return (np.fabs(mat1 - mat2)/mat1).mean()

def mat_MAE(mat1_file, mat2_file) -> float:
    mat1 = loadmat(mat1_file)
    mat2 = loadmat(mat2_file)
    return np.sum(np.abs(mat1 - mat2)) / np.prod(mat1.shape)

def mat_MaxErr(mat1_file, mat2_file) -> float:
    mat1 = loadmat(mat1_file)
    mat2 = loadmat(mat2_file)
    return np.max(np.abs(mat1 - mat2))

def list_files(directory):
    files = []
    try:
        with os.scandir(directory) as entries:
            for entry in entries:
                if entry.is_file():
                    files.append(entry.name)
    except FileNotFoundError:
        print(f"The directory '{directory}' does not exist.")
    except PermissionError:
        print(f"Permission denied to access the directory '{directory}'.")

    files.sort(key=lambda s: (len(s), s))
    return files

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <mat_dir1> <mat_dir2>")
        exit(1)

    dirname1 = sys.argv[1]
    dirname2 = sys.argv[2]
    dir1 = list_files(dirname1)
    dir2 = list_files(dirname2)

    dir1.reverse()
    dir2.reverse()

    while len(dir1) > 0 or len(dir2) > 0:
        file1 = dir1.pop()
        file2 = dir2.pop()

        if file1 != file2:
            print(f"'{dirname1}/{file1}' does not match '{dirname2}/{file2}'")
            continue

        file1 = f"{dirname1}/{file1}"
        file2 = f"{dirname2}/{file2}"

        mat1 = loadmat(file1)
        mat2 = loadmat(file2)

        print(f"====={file1} vs {file2}=====")
        print(f"Max Err: {mat_MaxErr(file1, file2)}")
        print(f"MAE: {mat_MAE(file1, file2)}")
        print(f"MAPE: {mat_MAPE(file1, file2)}")