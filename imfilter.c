
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <netpbm/pam.h>

typedef struct pixel_struct {
  int r;
  int g;
  int b;
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
    printf("reading row %d\n", i+1);
    pnm_readpamrow(pamImage, tuplerow);
    printf("read %d row\n", i+1);
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

pixel_t ** convolve(pixel_t **image, struct pam * pamImage, double **kernel, struct pam * pamKernel) {
  int i, j, k, ik, jk;
  int center_x = pamKernel->width/2;
  int center_y = pamKernel->height/2;
  double temp_result;

  pixel_t ** result = (pixel_t**) malloc(sizeof(pixel_t*)*(pamImage->height));
  for(i = center_y; i < pamImage->height-center_y; i++) {
    result[i] = (pixel_t*) malloc(sizeof(pixel_t)*width);
    for (j = center_x; j < pamImage->width-center_x; j++) {
      temp_result = 0;
      for (ik = 0; ik < pamKernel->height; ik++) {
        for (jk = 0; jk < pamKernel->width; jk++) {
          temp_result += image[i-center_y+ik][j - center_x + jk]*kernel[ik,jk];
        }
      }
      result[i][j] = temp_result;
    }
  }
  return result;
}

int main(int argc, char **argv) {
  struct pam pamImage, pamMask, pamOutput;
  pm_init(argv[0],0);
  pixel_t ** imArray;
  FILE * imageFile = fopen("StopSign2.ppm", "r");
  FILE * maskFile = fopen("gaussian.pgm","r");
  //FILE * output = fopen("output.ppm", "w");
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  pnm_readpaminit(maskFile, &pamMask, PAM_STRUCT_SIZE(tuple_type));
  //pamOutput = pamImage;
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

  freeRGBImage(imArray, &pamImage);


  //free(imArray);

  fclose(imageFile);
  fclose(maskFile);

  return 0;
}
