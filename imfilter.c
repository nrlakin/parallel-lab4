#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
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

pixel_t ** loadRGBImage(struct pam * pamImage) {
  tuple *tuplerow;
  pixel_t temp_pixel;
  int i, j;
  int width = pamImage->width;
  int height = pamImage->height;
  int depth = pamImage->depth;
  tuplerow = pnm_allocpamrow(pamImage);
  printf("allocating image...\n");
  pixel_t ** im = (pixel_t**) malloc(sizeof(pixel_t*)*height);
  //pixel_t ** im = (pixel_t**) malloc(sizeof(pixel_t)*width*height);
  printf("one pixel is %d bytes\n", sizeof(pixel_t));
  printf("initializing image...\n");
  for (i = 0; i < height; i++) {
    pnm_readpamrow(pamImage, tuplerow);
    im[i] = (pixel_t*)malloc(sizeof(pixel_t)*width);
    for (j = 0; j < width; j++) {
      im[i][j].r = tuplerow[j][0];
      im[i][j].g = tuplerow[j][1];
      im[i][j].b = tuplerow[j][2];
    }
  }
  pnm_freepamrow(tuplerow);
  printf("allocated image.\n");
  return im;
}

void freeRGBImage(pixel_t ** image, struct pam * pamImage) {
  int i;
  for (i = 0; i < pamImage->height; i++) {
    free(image[i]);
  }
  free(image);
}

double ** loadKernel(struct pam * pamImage) {
  tuple *tuplerow;
  int i, j;
  int width = pamImage->width;
  int height = pamImage->height;
  int maxval = pamMask->maxval;
  tuplerow = pnm_allocpamrow(pamImage);
  double ** im = (double**) malloc(sizeof(double*)*height);
  for (i = 0; i < height; i++) {
    pnm_readpamrow(pamImage, tuplerow);
    im[i] = (double*)malloc(sizeof(double)*width);
    for (j = 0; j < width; j++) {
      im[i][j] = (double) ((8*tuplerow[j][0]/maxval) - 4);
    }
  }
  pnm_freepamrow(tuplerow);
  return im;
}

void freeKernel(double ** image, struct pam * pamImage) {
  int i;
  for (i = 0; i < pamImage->height; i++) {
    free(image[i]);
  }
  free(image);
}

pixel_t ** convolve(pixel_t **image, struct pam * pamImage, double **kernel, struct pam * pamKernel) {
  int i, j, k, ik, jk;
  int center_x = pamKernel->width/2;
  int center_y = pamKernel->height/2;
  double temp_result_r, temp_result_g, temp_result_b;

  pixel_t ** result = (pixel_t**) malloc(sizeof(pixel_t*)*(pamImage->height));
  for(i = center_y; i < pamImage->height-center_y; i++) {
    result[i] = (pixel_t*) malloc(sizeof(pixel_t)*pamImage->width);
    for (j = center_x; j < pamImage->width-center_x; j++) {
      temp_result_r = 0;
      temp_result_g = 0;
      temp_result_b = 0;
      for (ik = 0; ik < pamKernel->height; ik++) {
        for (jk = 0; jk < pamKernel->width; jk++) {
          temp_result_r += image[i-center_y+ik][j - center_x + jk].vector[0]*kernel[ik][jk];
          temp_result_g += image[i-center_y+ik][j - center_x + jk].vector[1]*kernel[ik][jk];
          temp_result_b += image[i-center_y+ik][j - center_x + jk].vector[2]*kernel[ik][jk];
        }
      }
      result[i][j].vector[0] = temp_result_r;
      result[i][j].vector[1] = temp_result_g;
      result[i][j].vector[2] = temp_result_b;
    }
  }
  return result;
}

void writeMatrixToFile(struct pam * outPam, pixel_t **image, struct pam * imgPam) {
  int i, j;
  tuple * tuplerow;
  pnm_writepaminit(outPam);
  tuplerow = pnm_allocpamrow(imgPam);
  for (i = 0; i < imgPam->height; i++) {
      pnm_readpamrow(imgPam, tuplerow);
      for (j = 0; j < imgPam->width; j++) {
        tuplerow[i][j][0] = pixel_t[i][j].r;
        tuplerow[i][j][1] = pixel_t[i][j].g;
        tuplerow[i][j][2] = pixel_t[i][j].b;
      }
      pnm_writepamrow(outPam, tuplerow);
  }
  pnm_freepamrow(tuplerow);
}

int main(int argc, char **argv) {
  struct pam pamImage, pamMask, pamOutput;
  pm_init(argv[0],0);
  pixel_t ** imArray;
  double ** kerArray;
  FILE * imageFile = fopen("StopSign2.ppm", "r");
  FILE * maskFile = fopen("gaussian.pgm","r");
  FILE * output = fopen("output.ppm", "w");
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  pnm_readpaminit(maskFile, &pamMask, PAM_STRUCT_SIZE(tuple_type));
  pamOutput = pamImage;
  pamOutput.file = output;
  //printf("format: %s\n", inpam.format);
  printf("size: %d\n", pamImage.size);
  printf("len: %d\n", pamImage.len);
  printf("height: %d\n", pamImage.height);
  printf("width: %d\n", pamImage.width);
  printf("depth: %d\n", pamImage.depth);
  printf("max val: %d\n", pamImage.maxval);

  printf("mask info: \n");
  printf("size: %d\n", pamMask.size);
  printf("len: %d\n", pamMask.len);
  printf("height: %d\n", pamMask.height);
  printf("width: %d\n", pamMask.width);
  printf("depth: %d\n", pamMask.depth);
  printf("max val: %d\n", pamMask.maxval);

  printf("row 1 values:\n");
  imArray = loadRGBImage(&pamImage);
  kerArray = loadKernel(&pamMask);

  freeRGBImage(imArray, &pamImage);
  freeKernel(kerArray, &pamMask);

  fclose(imageFile);
  fclose(maskFile);
  fclose(output);

  return 0;
}
