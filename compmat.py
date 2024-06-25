import sys, struct

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} <mat1> <mat2>")
    exit(1)

with open(sys.argv[1], 'rb') as f1:
    with open(sys.argv[2], 'rb') as f2:
        m1_x, m1_y, m1_z = struct.unpack('<QQQ', f1.read(24))
        m2_x, m2_y, m2_z = struct.unpack('<QQQ', f2.read(24))

        if (m1_x, m1_y, m1_z) != (m2_x, m2_y, m2_z):
            print(999999999.99)
            exit(1)

        maxerr = -1000000000.0
        
        for i in range(m1_x * m1_y * m1_z):
            val1, = struct.unpack('<d', f1.read(8))
            val2, = struct.unpack('<d', f2.read(8))

            maxerr = max(maxerr, abs(val1 - val2))

print(maxerr)
