
def getWindow(n_proc, proc_id, kern_width, im_width, im_height):
    n_rows = int(n_proc**.5)
    n_columns = n_proc/n_rows

    left_top_pad = getSmallPad(kern_width)
    right_bottom_pad = int(kern_width/2)
    x_offset = int(im_width/n_columns)
    y_offset = int(im_height/n_rows)
    row_index = int(proc_id/n_rows)
    col_index = proc_id%n_columns

    i_start = row_index * y_offset;
    if row_index == (n_rows-1):
        i_end = getTotal(kern_width, im_height)
    else:
        i_end = (row_index + 1) * y_offset + left_top_pad + right_bottom_pad;

    j_start = col_index * x_offset
    if col_index == (n_columns-1):
        j_end = getTotal(kern_width, im_width)
    else:
        j_end = (col_index + 1) * x_offset + left_top_pad + right_bottom_pad;

    return i_start, i_end, j_start, j_end

def getSmallPad(MaskDim):
  pad = int(MaskDim / 2) - (1-(MaskDim%2))
  if (pad < 0):
      pad = 0
  return pad

def getTotal(MaskDim, ImgDim):
  smallPad = getSmallPad(MaskDim)
  bigPad = MaskDim / 2
  return ImgDim + smallPad + bigPad
