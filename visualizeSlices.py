import os
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt

def loadmat(filename: str) -> np.array:
    with open(filename, 'rb') as f:
        dim_x, dim_y, dim_z = struct.unpack('<QQQ', f.read(24))
        data = np.frombuffer(f.read(), dtype=np.float64).reshape(dim_x, dim_y, dim_z)
        return data

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

    return files

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} <n_sections> <input_directory> <output_directory>")
        exit(1)

    n_sections = sys.argv[1]
    input_directory = sys.argv[2]
    output_directory = sys.argv[3]
    files = list_files(input_directory)
    files.sort(key=lambda s: (len(s), s))

    for file in files:
        f = f'{input_directory}/{file}'
        print(f)
        mat = loadmat(f)
        N1, N2, N3 = mat.shape

        X, Y, Z = np.meshgrid(np.arange(N1), np.arange(N2), -np.arange(N3))

        # figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        vmin = -2
        vmax = 2

        kw = {
            'vmin': vmin,
            'vmax': vmax,
            'levels': np.linspace(vmin, vmax, 100),
            'cmap': "coolwarm"
        }

        C = None
        for offset in np.linspace(0, N3-1, n_sections):
            offset = round(offset)
            C = ax.contourf(
                X[:, :, offset], Y[:, :, offset], mat[:, :, offset],
                zdir='z', offset=-offset, **kw
            )

        xmin, xmax = X.min(), X.max()
        ymin, ymax = Y.min(), Y.max()
        zmin, zmax = Z.min(), Z.max()
        ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])
        
        # account for LHS-RHS difference
        ax.set(
            xlabel='Y',
            ylabel='X',
            zlabel='Z'
        )

        ax.view_init(10, -30, 0)
        ax.set_box_aspect(None, zoom=0.9)

        fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Temperature')

        plt.savefig(f"{output_directory}/{file}.png")
        plt.close()