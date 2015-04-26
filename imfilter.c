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
  for (i = 0; i < height; i++) {
    printf("%.4f    %.4f    %.4f    %.4f    %.4f\n", im[i][0], im[i][1], im[i][2], im[i][3], im[i][4]);
  }
  return im;
}

void freeKernel(double ** image, struct pam * pamImage) {
  int i;
  for (i = 0; i < pamImage->height; i++) {
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

pixel_t ** convolve(pixel_t **image, struct pam * pamImage, double **kernel, struct pam * pamKernel) {
  int i, j, k, ik, jk, target_i, target_j;
  double ksum = 0;
  int center_x = pamKernel->width/2;
  int center_y = pamKernel->height/2;
  printf("kernel center_x: %d\n", center_x);
  printf("kernel center_y: %d\n", center_y);
  double temp_result_r, temp_result_g, temp_result_b;

  pixel_t ** result = (pixel_t**) malloc(sizeof(pixel_t*)*(pamImage->height));
  for(i = 0; i < pamImage->height; i++) {
    result[i] = (pixel_t*) malloc(sizeof(pixel_t)*pamImage->width);
    for (j = 0; j < pamImage->width; j++) {
      temp_result_r = 0;
      temp_result_g = 0;
      temp_result_b = 0;
      ksum = 0;
      for (ik = 0; ik < pamKernel->height; ik++) {
        for (jk = 0; jk < pamKernel->width; jk++) {
          target_i = i + (ik - center_y);
          target_j = j + (jk - center_x);
          if (target_i >= 0 && target_i < pamImage->height && target_j >= 0 && target_j < pamImage->width) {
            temp_result_r += ((double)image[target_i][target_j].vector[0])*kernel[ik][jk];
            temp_result_g += ((double)image[target_i][target_j].vector[1])*kernel[ik][jk];
            temp_result_b += ((double)image[target_i][target_j].vector[2])*kernel[ik][jk];
          }
          ksum += kernel[ik][jk];
        }
      }
      result[i][j].r = getPixel(temp_result_r, ksum, pamImage->maxval);
      result[i][j].g = getPixel(temp_result_g, ksum, pamImage->maxval);
      result[i][j].b = getPixel(temp_result_b, ksum, pamImage->maxval);
    }
  }
  return result;
}

void writeMatrixToFile(struct pam * outPam, pixel_t **image, struct pam * imgPam) {
  int i, j;
  tuple * tuplerow;
  pnm_writepaminit(outPam);
  tuplerow = pnm_allocpamrow(imgPam);
  printf("Writing output to file.\n");
  for (i = 0; i < imgPam->height; i++) {
      //printf("reading row %d\n", i);
      pnm_readpamrow(imgPam, tuplerow);
      //printf("read row %d\n", i);
      for (j = 0; j < imgPam->width; j++) {
        tuplerow[j][0] = image[i][j].r;
        tuplerow[j][1] = image[i][j].g;
        tuplerow[j][2] = image[i][j].b;
      }
      pnm_writepamrow(outPam, tuplerow);
  }
  pnm_freepamrow(tuplerow);
}

int main(int argc, char **argv) {
  struct pam pamImage, pamMask, pamOutput;
  pm_init(argv[0],0);
  pixel_t ** imArray, ** outputImage;
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

  imArray = loadRGBImage(&pamImage);
  kerArray = loadKernel(&pamMask);

  outputImage = convolve(imArray, &pamImage, kerArray, &pamMask);
  // CLUNKY
  rewind(imageFile);
  pnm_readpaminit(imageFile, &pamImage, PAM_STRUCT_SIZE(tuple_type));
  writeMatrixToFile(&pamOutput, outputImage, &pamImage);
  //writeMatrixToFile(&pamOutput, outputImage, &pamImage);


  freeRGBImage(imArray, &pamImage);
  freeRGBImage(outputImage, &pamOutput);
  freeKernel(kerArray, &pamMask);

  fclose(imageFile);
  fclose(maskFile);
  fclose(output);

  return 0;
}
