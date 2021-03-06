#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <netpbm/pam.h>

/*** Structure representing a single pixel. Values are 0-255; use ints to
  prevent overflow during normalization.  ***/
typedef struct pixel_struct {
  union {
    struct {
        int r;
        int g;
        int b;
    };
    int vector[3];
  };
} pixel_t;

/*** Structure for keeping track of image subwindow geometry. ***/
typedef struct {
  int i_start;
  int i_end;
  int j_start;
  int j_end;
  int l_pad;
  int r_pad;
  int t_pad;
  int b_pad;
} window_t;

/***
  void getWindow(window_t * destPtr, int n_proc, int proc_id, struct pam * pamKernel, struct pam * pamImage)

  Populate window_t structure *destPtr with window parameters given processor
  id, number of processors, and image/kernel geometries.
***/
void getWindow(window_t * destPtr, int n_proc, int proc_id, struct pam * pamKernel, struct pam * pamImage) {
    int n_rows = (int)sqrt(n_proc);
    while(n_proc%n_rows)n_rows--;
    int n_columns = n_proc/n_rows;
    int kern_width = pamKernel->width;
    int kern_height = pamKernel->height;
    int i_start, i_end, j_start, j_end;
    int left_pad = getSmallPad(kern_width);
    int top_pad = getSmallPad(kern_height);
    int right_pad = kern_width/2;
    int bottom_pad = kern_height/2;
    int x_offset = pamImage->width/n_columns;
    int y_offset = pamImage->height/n_rows;
    int row_index = proc_id/(n_proc/n_rows);
    int col_index = proc_id%n_columns;

    i_start = row_index * y_offset;
    i_end = (row_index + 1) * y_offset + top_pad + bottom_pad;
    i_end = (row_index == n_rows-1) ? getTotal(kern_height, pamImage->height) : i_end;
    destPtr->i_start = i_start;
    destPtr->i_end = i_end;

    j_start = col_index * x_offset;
    j_end = (col_index + 1) * x_offset + left_pad + right_pad;
    j_end = (col_index == n_columns - 1) ? getTotal(kern_width, pamImage->width) : j_end;
    destPtr->j_start = j_start;
    destPtr->j_end = j_end;
    destPtr->l_pad = left_pad;
    destPtr->r_pad = right_pad;
    destPtr->t_pad = top_pad;
    destPtr->b_pad = bottom_pad;
}

/***
  int getSmallPad(int MaskDim)

  Helper function to get right/bottom pad given x/y dimension.

***/
int getSmallPad(int MaskDim) {
  int pad = (MaskDim / 2) - (1-(MaskDim%2));
  if (pad < 0) pad = 0;
  return pad;
}

/***
  int getTotal(int MaskDim, int ImgDim)

  Helper function to calculate total dimensions of image after adding '0'
  padding.

***/
int getTotal(int MaskDim, int ImgDim) {
  int smallPad = getSmallPad(MaskDim);
  int bigPad = MaskDim / 2;
  return ImgDim + smallPad + bigPad;
}

/***
  pixel_t ** loadRGBImage(struct pam * pamImage, struct pam * pamMask)

  Return 2-D array of pixels, with a padding of zeroes around the image, given
  pam structures for the image and kernel.
***/
pixel_t ** loadRGBImage(struct pam * pamImage, struct pam * pamMask) {
  tuple *tuplerow;
  pixel_t temp_pixel;
  int i, j;
  int leftPad = getSmallPad(pamMask->width);
  int topPad = getSmallPad(pamMask->height);
  int width = getTotal(pamMask->width, pamImage->width);
  int height = getTotal(pamMask->height, pamImage->height);

  tuplerow = pnm_allocpamrow(pamImage);
  // printf("allocating image...\n");
  pixel_t ** im = (pixel_t**) malloc(sizeof(pixel_t*)*height);
  // printf("initializing image...\n");
  for (i = 0; i < height; i++) {
    im[i] = (pixel_t*)malloc(sizeof(pixel_t)*width);
    if ((topPad <= i) && (i < (topPad + pamImage->height))) {
      pnm_readpamrow(pamImage, tuplerow);
    }
    for (j = 0; j < width; j++) {
      if ((topPad <= i) && (i < (topPad + pamImage->height)) && (leftPad <= j) && (j < (leftPad + pamImage->width))) {
        im[i][j].r = tuplerow[j-leftPad][0];
        im[i][j].g = tuplerow[j-leftPad][1];
        im[i][j].b = tuplerow[j-leftPad][2];
      } else {
        im[i][j].r = 0;
        im[i][j].g = 0;
        im[i][j].b = 0;
      }
    }
  }
  pnm_freepamrow(tuplerow);
  // printf("allocated image.\n");
  return im;
}

/***
  double ** loadKernel(struct pam * pamImage)

  Return a 2D array of doubles representing the convolution kernel, given a
  pam struct representing the pgm file.
***/
double ** loadKernel(struct pam * pamImage) {
  tuple *tuplerow;
  int i, j;
  int width = pamImage->width;
  int height = pamImage->height;
  int maxval = pamImage->maxval;
  tuplerow = pnm_allocpamrow(pamImage);

  double ** im = (double**) malloc(sizeof(double*)*height);
  for (i = 0; i < height; i++) {
    pnm_readpamrow(pamImage, tuplerow);
    im[i] = (double*)malloc(sizeof(double)*width);
    for (j = 0; j < width; j++) {
      im[i][j] = ((8*((double)tuplerow[j][0])/maxval) - 4);
    }
  }
  pnm_freepamrow(tuplerow);
  //for (i = 0; i < height; i++) {
  //  printf("%.4f    %.4f    %.4f    %.4f    %.4f\n", im[i][0], im[i][1], im[i][2], im[i][3], im[i][4]);
  //}
  return im;
}

/***
  void freeRGBImage(pixel_t ** image, int height)

  Helper function to free memory for 2D array representing image.
***/
void freeRGBImage(pixel_t ** image, int height) {
  int i;
  for (i = 0; i < height; i++) {
    free(image[i]);
  }
  free(image);
}

/***
  void freeKernel(double ** image, int height)

  Helper function to free memory for 2D array representing kernel.
***/
void freeKernel(double ** image, int height) {
  int i;
  for (i = 0; i < height; i++) {
    free(image[i]);
  }
  free(image);
}

/***
  int getPixel(double temp_pixel, double ksum, int maxval)

  Normalize (roughly) a single pixel after a convolution.
***/
int getPixel(double temp_pixel, double ksum, int maxval) {
  if (ksum <= 0) ksum = 1;
  temp_pixel /= ksum;
  if (temp_pixel > maxval) {
    temp_pixel = maxval;
  } else if (temp_pixel < 0){
    temp_pixel = 0;
  }
  return (int) temp_pixel;
}

/***
  pixel_t ** convolve(pixel_t **image, window_t * window, int img_maxval, double **kernel, struct pam * pamKernel)

  Given 2D arrays for image (including padding) and kernel, return a pointer to
  a 2D array representing the output.  Input can be a whole image or a subwindow
  of the image.  Geometric data about the image is passed in window_t; geometric
  information about the kernel is included in the pamKernel structure.  Note
  that the output image is not padded in any way and so may have different
  dimensions than the input image.
***/
pixel_t ** convolve(pixel_t **image, window_t * window, int img_maxval, double **kernel, struct pam * pamKernel) {
  int i, j, k, ik, jk, target_i, target_j;
  int width = window->j_end - window->j_start;
  int height = window->i_end - window->i_start;
  int img_height = height - (window->t_pad + window->b_pad);
  int img_width = width - (window->l_pad + window->r_pad);
  double ksum = 0;
  int center_x = (pamKernel->width-1)/2;
  int center_y = (pamKernel->height-1)/2;
  // printf("kernel center_x: %d\n", center_x);
  // printf("kernel center_y: %d\n", center_y);
  double temp_result_r, temp_result_g, temp_result_b;

  //int width = getTotal(pamKernel->width, img_width);
  //int height = getTotal(pamKernel->height, img_height);
  int leftPad = window->l_pad;
  int topPad = window->t_pad;

  // if(omp_get_thread_num()==0) {
  //   printf("leftPad=%d\n",leftPad);
  // }
  pixel_t ** result = (pixel_t**) malloc(sizeof(pixel_t*)*img_height);
  for(i = topPad; i < (topPad + img_height); i++) {
    result[i-topPad] = (pixel_t*) malloc(sizeof(pixel_t)*img_width);
    for (j = leftPad; j < (leftPad + img_width); j++) {
      //printf("working on px: %d, %d\n", i, j);
      temp_result_r = 0;
      temp_result_g = 0;
      temp_result_b = 0;
      ksum = 0;
      for (ik = 0; ik < pamKernel->height; ik++) {
        for (jk = 0; jk < pamKernel->width; jk++) {
          target_i = i + (ik - center_y);
          target_j = j + (jk - center_x);
          if (target_i >= 0 && target_i < height && target_j >= 0 && target_j < width) {
            temp_result_r += ((double)image[target_i][target_j].vector[0])*kernel[ik][jk];
            temp_result_g += ((double)image[target_i][target_j].vector[1])*kernel[ik][jk];
            temp_result_b += ((double)image[target_i][target_j].vector[2])*kernel[ik][jk];
          }
          ksum += kernel[ik][jk];
        }
      }
      result[i-topPad][j-leftPad].r = getPixel(temp_result_r, ksum, img_maxval);
      result[i-topPad][j-leftPad].g = getPixel(temp_result_g, ksum, img_maxval);
      result[i-topPad][j-leftPad].b = getPixel(temp_result_b, ksum, img_maxval);
    }
  }
  return result;
}

/***
  void writeMatrixToFile(struct pam * outPam, pixel_t **image, struct pam * imgPam, int height, int width)

  Given pointer to output pam structure (which implies an output ppm file) and
  a 2D image array, write output to file.
***/
void writeMatrixToFile(struct pam * outPam, pixel_t **image, struct pam * imgPam, int height, int width) {
  int i, j;
  tuple * tuplerow;
  pnm_writepaminit(outPam);
  tuplerow = pnm_allocpamrow(imgPam);
  // printf("Writing output to file.\n");
  for (i = 0; i < height; i++) {
      //printf("reading row %d\n", i);
      pnm_readpamrow(imgPam, tuplerow);
      //printf("read row %d\n", i);
      for (j = 0; j < width; j++) {
        tuplerow[j][0] = image[i][j].r;
        tuplerow[j][1] = image[i][j].g;
        tuplerow[j][2] = image[i][j].b;
      }
      pnm_writepamrow(outPam, tuplerow);
  }
  pnm_freepamrow(tuplerow);
}

/***
  pixel_t ** copyImgMatrix(pixel_t ** whole_img, window_t * window)

  Copy a subwindow of an image to a new 2D array.  Rather than have all
  processors read the same global array, copy subwindows to local arrays and
  write output to shared memory after convolving.
***/
pixel_t ** copyImgMatrix(pixel_t ** whole_img, window_t * window) {
  int i, j;
  // printf("window gotten: %d, %d, %d, %d\n", window->i_start, window->i_end,window->j_start,window->j_end);
  int height = (window->i_end - window->i_start);
  int width = (window->j_end - window->j_start);
  // printf("mallocing window rows\n");
  pixel_t ** sub_img = (pixel_t**) malloc(sizeof(pixel_t*)*(height));
  //printf ("malloc'd rows on %d\n", omp_get_thread_num());
  for (i = window->i_start; i < window->i_end; i++) {
    sub_img[i-window->i_start] = (pixel_t*)malloc(sizeof(pixel_t)*width);
    for (j = window->j_start; j < window->j_end; j++) {
        sub_img[i-window->i_start][j-window->j_start].r = whole_img[i][j].r;
        sub_img[i-window->i_start][j-window->j_start].g = whole_img[i][j].g;
        sub_img[i-window->i_start][j-window->j_start].b = whole_img[i][j].b;
    }
  }
  return sub_img;
}

/***
  void copyResult(pixel_t ** subImage, pixel_t ** destImage, window_t * window)

  Helper function to copy result subwindow back to shared memory.
***/
void copyResult(pixel_t ** subImage, pixel_t ** destImage, window_t * window) {
  int i_start = window->i_start + window->t_pad;
  int i_end = window->i_end - window->b_pad;
  int j_start = window->j_start + window->l_pad;
  int j_end = window->j_end - window->r_pad;
  int i, j, k;
  for (i = i_start; i < i_end; i++) {
    for (j = j_start; j < j_end; j++) {
      //for (k = 0; k < 3; k++) {
        //destImage[i][j].vector[k] = subImage[i-i_start][j-j_start].vector[k];
        destImage[i][j].r = subImage[i-i_start][j-j_start].r;
        destImage[i][j].g = subImage[i-i_start][j-j_start].g;
        destImage[i][j].b = subImage[i-i_start][j-j_start].b;
      //}
    }
  }
}

/*** Main function. ***/
int main(int argc, char **argv) {
  int opt;
  struct pam pamImage, pamMask, pamOutput;
  pm_init(argv[0],0);
  pixel_t ** imArray, ** outputImage;
  double ** kerArray;
  char *infile, *outfile;

  // Very brittle reading in command line args...
  FILE * imageFile = fopen(argv[1], "r");
  FILE * output = fopen(argv[2], "w");
  FILE * maskFile = fopen(argv[3],"r");
  int n_threads = atoi(argv[4]);
  int n_iter = atoi(argv[5]);
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  pnm_readpaminit(maskFile, &pamMask, PAM_STRUCT_SIZE(tuple_type));
  pamOutput = pamImage;
  pamOutput.file = output;

  kerArray = loadKernel(&pamMask);
  imArray = loadRGBImage(&pamImage, &pamMask);
  // printf("finished loading\n");


  // Begin parallel part.
  #pragma omp parallel num_threads(n_threads) shared(imArray, n_iter)
  {
    window_t my_window;
    pixel_t ** sub_im, ** sub_out;
    int height, width, totalpad, t;

    int my_thread = omp_get_thread_num();
    int n_proc = omp_get_num_threads();
    // printf("I am thread %d of %d\n", my_thread, n_proc);
    getWindow(&my_window, n_threads, my_thread, &pamMask, &pamImage);

    // Make sure everyone reads before writing output...
    #pragma omp barrier
    //printf("thread %d got window.\n", my_thread);
    for(t = 0; t < n_iter; t++) {
      sub_im = copyImgMatrix(imArray, &my_window);
      //printf("window copied\n");
      height = (my_window.i_end - my_window.i_start);
      width = (my_window.j_end - my_window.j_start);
      sub_out = convolve(sub_im, &my_window, pamImage.maxval, kerArray, &pamMask);
      printf("got through convolution...\n");

    // Only one processor writes to common output at a time
    #pragma omp critical
    {
      copyResult(sub_out, imArray, &my_window);
    }

      freeRGBImage(sub_im, height);
      freeRGBImage(sub_out, height-(my_window.t_pad+my_window.b_pad));

    // Let everyone finish before next iteration. Not pretty.
    #pragma omp barrier
    }
  }

  // CLUNKY
  rewind(imageFile);
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  writeMatrixToFile(&pamOutput, imArray, &pamImage, pamImage.height, pamImage.width);

  freeRGBImage(imArray, getTotal(pamMask.height, pamImage.height));
  // freeRGBImage(outputImage, pamImage.height);
  freeKernel(kerArray, pamMask.height);

  fclose(imageFile);
  fclose(maskFile);
  fclose(output);

  return 0;
}
