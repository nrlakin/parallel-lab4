import numpy as np


def kernelToPGM(arr, filename='kernel.pgm'):
    width = arr.shape[1]
    height = arr.shape[0]
    max_val = 255
    header = bytearray('P5\x0A'+str(width)+' '+str(height)+'\x0A' + str(max_val) + '\x0A')
    pix_map = ((arr.astype(np.float64) + 4)*max_val/8).astype(np.uint8)
    with open(filename, 'w') as f:
        f.write(header)
        for i in range(0, pix_map.shape[0]):
            f.write(bytearray(list(pix_map[i, :])))


def PGMtoKernel(filename='kernel.pgm'):
    with open(filename, 'r') as f:
        contents = f.read()
    head = contents.split()
    width = int(head[1])
    height = int(head[2])
    max_val = int(head[3])
    n_pix = width*height
    pix_map = np.array(bytearray(contents[-n_pix:])).reshape(height, width).astype(np.float64)
    result = 8*pix_map/max_val - 4
    return result
