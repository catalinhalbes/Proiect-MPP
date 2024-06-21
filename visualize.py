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
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <input_directory> <output_directory>")
        exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    files = list_files(input_directory)
    files.sort(key=lambda s: (len(s), s))

    for file in files:
        f = f'{input_directory}/{file}'
        print(f)
        mat = loadmat(f)
        N1, N2, N3 = mat.shape

        X, Y, Z = np.meshgrid(np.arange(N1), np.arange(N2), np.arange(N3))

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

        _ = ax.contourf(
            X[:, :, 0], Y[:, :, 0], mat[:, :, -1],
            zdir='z', offset=Z.max(), **kw
        )
        _ = ax.contourf(
            X[0, :, :], mat[0, :, :], Z[0, :, :],
            zdir='y', offset=0, **kw
        )
        C = ax.contourf(
            mat[:, -1, :], Y[:, -1, :], Z[:, -1, :],
            zdir='x', offset=X.max(), **kw
        )

        xmin, xmax = X.min(), X.max()
        ymin, ymax = Y.min(), Y.max()
        zmin, zmax = Z.min(), Z.max()
        ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

        edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
        ax.plot([xmax, xmax], [ymin, ymax], [zmax, zmax], **edges_kw)
        ax.plot([xmin, xmax], [ymin, ymin], [zmax, zmax], **edges_kw)
        ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)
        
        ax.set(
            xlabel='Y',
            ylabel='X',
            zlabel='Z'
        )

        ax.view_init(40, -30, 0)
        ax.set_box_aspect(None, zoom=0.9)

        fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Temperature')

        plt.savefig(f"{output_directory}/{file}.png")
        plt.close()