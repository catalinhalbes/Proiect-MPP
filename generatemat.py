import sys, struct

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print(f"Usage: python {sys.argv[0]} <x dim> <y dim> <z dim> <left/hot side temp> <middle temp> <right/cold side temp> <out_file>")
        exit(1)

    x_dim       = int(sys.argv[1])
    y_dim       = int(sys.argv[2])
    z_dim       = int(sys.argv[3])
    left_temp   = float(sys.argv[4])
    mid_temp    = float(sys.argv[5])
    right_temp  = float(sys.argv[6])
    out_file    = sys.argv[7]

    with open(out_file, 'wb') as f:
        f.write(struct.pack('<Q', x_dim))
        f.write(struct.pack('<Q', y_dim))
        f.write(struct.pack('<Q', z_dim))
        
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if i == 0:
                        f.write(struct.pack('<d', left_temp))
                    elif i == x_dim - 1:
                        f.write(struct.pack('<d', right_temp))
                    else:
                        f.write(struct.pack('<d', mid_temp))
                    