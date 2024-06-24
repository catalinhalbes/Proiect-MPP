import sys, struct

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: python {sys.argv[0]} <x dim> <y dim> <z dim> <out_file>")
        exit(1)

    x_dim       = int(sys.argv[1])
    y_dim       = int(sys.argv[2])
    z_dim       = int(sys.argv[3])
    out_file    = sys.argv[4]

    with open(out_file, 'wb') as f:
        f.write(struct.pack('<Q', x_dim))
        f.write(struct.pack('<Q', y_dim))
        f.write(struct.pack('<Q', z_dim))
        
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    f.write(struct.pack('<d', i * (y_dim * z_dim) + j * z_dim + k))
                    