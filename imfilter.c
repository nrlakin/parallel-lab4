#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <netpbm/pam.h>

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

typedef struct {
  int i_start;
  int i_end;
  int j_start;
  int j_end;
} window_t;

void getWindow(window_t * destPtr, int n_proc, int proc_id, struct pam * pamKernel, struct pam * pamImage) {
    int n_rows = (int)sqrt(n_proc);
    int n_columns = n_proc/n_rows;
    int kern_width = pamKernel->width;
    int kern_height = pamKernel->height;
    int i_start, i_end, j_start, j_end;
    int left_top_pad = getSmallPad(kern_width);
    int right_bottom_pad = kern_width/2;
    int x_offset = pamImage->width/n_columns;
    int y_offset = pamImage->height/n_rows;
    int row_index = (n_rows==1) ? 0 : proc_id/n_rows;
    int col_index = (n_columns == 1) ? 0 : proc_id%n_columns;

    i_start = row_index * y_offset;
    i_end = (row_index + 1) * y_offset + left_top_pad + right_bottom_pad;
    i_end = (row_index == n_rows-1) ? getTotal(kern_height, pamImage->height) : i_end;
    destPtr->i_start = i_start;
    destPtr->i_end = i_end;

    j_start = col_index * x_offset;
    j_end = (col_index + 1) * x_offset + left_top_pad + right_bottom_pad;
    j_end = (col_index == n_columns - 1) ? getTotal(kern_width, pamImage->width) : j_end;
    destPtr->j_start = j_start;
    destPtr->j_end = j_end;
}

int getSmallPad(int MaskDim) {
  int pad = (MaskDim / 2) - (1-(MaskDim%2));
  if (pad < 0) pad = 0;
  return pad;
}

int getTotal(int MaskDim, int ImgDim) {
  int smallPad = getSmallPad(MaskDim);
  int bigPad = MaskDim / 2;
  return ImgDim + smallPad + bigPad;
}

pixel_t ** loadRGBImage(struct pam * pamImage, struct pam * pamMask) {
  tuple *tuplerow;
  pixel_t temp_pixel;
  int i, j;
  int leftPad = getSmallPad(pamMask->width);
  int topPad = getSmallPad(pamMask->height);
  int width = getTotal(pamMask->width, pamImage->width);
  int height = getTotal(pamMask->height, pamImage->height);

  tuplerow = pnm_allocpamrow(pamImage);
  printf("allocating image...\n");
  pixel_t ** im = (pixel_t**) malloc(sizeof(pixel_t*)*height);
  printf("initializing image...\n");
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
  printf("allocated image.\n");
  return im;
}

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

void freeRGBImage(pixel_t ** image, int height) {
  int i;
  for (i = 0; i < height; i++) {
    free(image[i]);
  }
  free(image);
}

void freeKernel(double ** image, int height) {
  int i;
  for (i = 0; i < height; i++) {
    free(image[i]);
  }
  free(image);
}

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

pixel_t ** convolve(pixel_t **image, int img_height, int img_width, int img_maxval, double **kernel, struct pam * pamKernel) {
  int i, j, k, ik, jk, target_i, target_j;
  double ksum = 0;
  int center_x = pamKernel->width/2;
  int center_y = pamKernel->height/2;
  printf("kernel center_x: %d\n", center_x);
  printf("kernel center_y: %d\n", center_y);
  double temp_result_r, temp_result_g, temp_result_b;

  int width = getTotal(pamKernel->width, img_width);
  int height = getTotal(pamKernel->height, img_height);
  int leftPad = getSmallPad(pamKernel->width);
  int topPad = getSmallPad(pamKernel->height);

  pixel_t ** result = (pixel_t**) malloc(sizeof(pixel_t*)*img_height);
  for(i = topPad; i < (topPad + img_height); i++) {
    result[i-topPad] = (pixel_t*) malloc(sizeof(pixel_t)*img_width);
    for (j = leftPad; j < (leftPad + img_width); j++) {
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

void writeMatrixToFile(struct pam * outPam, pixel_t **image, struct pam * imgPam, int height, int width) {
  int i, j;
  tuple * tuplerow;
  pnm_writepaminit(outPam);
  tuplerow = pnm_allocpamrow(imgPam);
  printf("Writing output to file.\n");
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

pixel_t ** copyImgMatrix(pixel_t ** whole_img, window_t * window) {
  int i, j;
  printf("window gotten: %d, %d, %d, %d\n", window->i_start, window->i_end,window->j_start,window->j_end);
  int height = (window->i_end - window->i_start - 1);
  int width = (window->j_end - window->j_start - 1);
  printf("mallocing window rows\n");
  pixel_t ** sub_img = (pixel_t**) malloc(sizeof(pixel_t*)*(height));
  for (i = window->i_start; i < window->i_end; i++) {
    sub_img[i] = (pixel_t*)malloc(sizeof(pixel_t)*width);
    for (j = window->j_start; j < window->j_end; j++) {
        sub_img[i][j].r = whole_img[i][j].r;
        sub_img[i][j].g = whole_img[i][j].g;
        sub_img[i][j].b = whole_img[i][j].b;
    }
  }
  return sub_img;
}

int main(int argc, char **argv) {
  int opt;
  struct pam pamImage, pamMask, pamOutput;
  pm_init(argv[0],0);
  pixel_t ** imArray, ** outputImage;
  double ** kerArray;
  char *infile, *outfile;

  //while ((opt = getopt(argc, argv, "oin")) != -1) {
  //  switch(opt) {
  //    case 'o': outfile =
  //  }
  //}
  //FILE * imageFile = fopen("StopSign2.ppm", "r");
  FILE * imageFile = fopen(argv[1], "r");
  FILE * maskFile = fopen("gaussian.pgm","r");
  FILE * output = fopen(argv[2], "w");
  int n_threads = atoi(argv[3]);
  int n_iter = atoi(argv[4]);
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  pnm_readpaminit(maskFile, &pamMask, PAM_STRUCT_SIZE(tuple_type));
  pamOutput = pamImage;
  pamOutput.file = output;

  kerArray = loadKernel(&pamMask);
  imArray = loadRGBImage(&pamImage, &pamMask);
  printf("finished loading\n");

  window_t my_window;
  pixel_t ** sub_im;
  int height, width;

  #pragma omp parallel num_threads(n_threads) private(sub_im, my_window, height, width)
  {
    int my_thread = omp_get_thread_num();
    int n_proc = omp_get_num_threads();
    printf("I am thread %d of %d\n", my_thread, n_proc);
    getWindow(&my_window, n_threads, my_thread, &pamMask, &pamImage);
    sub_im = copyImgMatrix(imArray, &my_window);
    int height = (my_window.i_end - my_window.i_start - 1);
    int width = (my_window.j_end - my_window.j_start - 1);
    outputImage = convolve(sub_im, height, width, pamImage.maxval, kerArray, &pamMask);
    freeRGBImage(sub_im, height);
    printf("finished convolution\n");
  }
  // CLUNKY
  rewind(imageFile);
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  // writeMatrixToFile(&pamOutput, outputImage, &pamImage, pamImage.height, pamImage.width);

  freeRGBImage(imArray, getTotal(pamMask.height, pamImage.height));
  // freeRGBImage(outputImage, pamImage.height);
  freeKernel(kerArray, pamMask.height);

  fclose(imageFile);
  fclose(maskFile);
  fclose(output);

  return 0;
}
