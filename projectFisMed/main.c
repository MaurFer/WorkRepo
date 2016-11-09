/* Programa que lee una imagen en formato .pgm (ASCII) y  escribe otra ya procesada */
#define usageQT

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#define Maxline 1000
#define sizeLevel 256
#define PI 3.1415926536
#define dim3 3
#define dim5 5
#define dim7 7


void read_pgm_file(char *filename);
void write_pgm_file(char *filename);
void write_pgm_hist();
void process_pgm_file();
void WriteTextPoint(void);

float inter_lin(float,float);
float inter_Neighbor(float x,float y);
float inter_bicubic(float x,float y); //*

void getHistogram(int *inputPointer);
double generatorNoiseGauss(int m, float d);
double generatoNoisedoubleRandn (double mu, double sigma);
void equalizeImage();
void convolucionVector();

int   filterTransform(int valPixel);
float filterLog(int valPixel);
float filterExp(int valPixel);
void  filterGradient(int dim); //*
void  filterGaussian(int dim);  //*
void  filterLaplacian(int dim); //*


/* Variables globales */
int width, height, bit_depth;
int width_out, height_out, bit_depth_out;
int width_hist, height_hist, bit_depth_hist;
int *image_pointer, *image_temp, *image_histogram;
int *image_pointer_out;

float *sumHist;

float histogram[sizeLevel];
float histogramNormal[sizeLevel];

float eqMax;
float eqMin;


int main(int argc, char *argv[])
{
if(argc < 3) {
    //printf("error en los datos de entrada\n"); abort();
    printf("********************************************************\n");
    printf("Usage: imageA_input imageZ_Output operation param\n");
    printf("\n");
    printf("operation: 0  Graph histogram \n");
    printf("           1  Equalize image    | Param: MinLevel MaxLevel\n");
    printf("           2  Subtrac A-B       | param: imageB_input\n");
    printf("           3  Filter Binary:    | Param: Value \n");
    printf("           4  Filter Log        | Param: Value c \n");
    printf("           5  Filter Exp        | Param: Value c \n");
    printf("           6  Filter Gradiente  | Param: 3x3 5x5 7x7 \n");
    printf("           7  Filter Laplacian  | Param: 3x3 5x5 7x7 \n");
    printf("           8  Filter Gaussiano  | Param: 3x3 5x5 7x7 \n");
    printf("           9  Resize lineal     | Param: PixelWxH  \n");
    printf("           10 Resize bilineal   | Param: PixelWxH  \n");
    printf("           11 Resize bicubica   | Param: PixelWxH  \n");
    printf("           12 Add Noise         | \n");
    printf("********************************************************\n");
    printf("\n");
    printf("\n");
    //abort();

    argv[1] = "output_ImagenA";
    argv[2] = "ImagenZ.pgm";
    argv[3] = "7";
    argv[4] = "3";
    argv[5] = "1";
    argc = 4;
    eqMax = atoi(argv[4]);
    eqMin = atoi(argv[5]);
}

read_pgm_file(argv[1]);

int i;
char file[25];
char command[25];
//select operation

printf("%i", atoi(argv[3]));
switch (atoi(argv[3])){

//This operation Equalize image
case 1:
    if(argc < 3){
        printf("##################################################\n");
        printf("operation Equalize | ERR : input Image\n");
        printf("\n");
        abort();
    }

    strcpy(file, "output_");
    strncat(file, argv[1], strlen(argv[1])-4);
    printf("\nequalize %s to %s: \n\n", argv[1], file);

    //initialize
    width_out=width;
    height_out=height;
    bit_depth_out=bit_depth;

    width_hist = sizeLevel;
    height_hist = sizeLevel;
    bit_depth_hist=1;
    image_histogram=(int*)malloc(sizeof(int)*width_hist*height_hist);
    sumHist=(int*)malloc(sizeof(int)*sizeLevel);

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

    getHistogram(image_pointer);
    equalizeImage();

    write_pgm_file(file);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, file, strlen(file));
    strncat(command, " \n", 2);
    break;

//This operation substration image A - image B
case 2:
    if(argc < 4){
        printf("ERR: input ImageA - input ImageB\n");
        abort();
    }


    printf("\nOperation subtrac  %s - %s: \n\n", argv[1], argv[4]);

    //initialize
    width_out=width;
    height_out=height;
    bit_depth_out=bit_depth;
    image_temp=(int *)malloc(sizeof(int)*width_out*height_out);

    for(i=0; i< width_out* height_out; i++)
        image_temp[i] = image_pointer[i];  //imageA (argv[1])

    read_pgm_file(argv[4]);     //imageA (argv[4])

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

    int pixel_value;
    for(i=0; i< width_out* height_out; i++){
        pixel_value =image_temp[i]- image_pointer[i];
        //if(pixel_value <0) pixel_value = pixel_value*(-1);
        if(pixel_value <0) pixel_value =0;
        image_pointer_out[i] = pixel_value;
    }


    write_pgm_file(argv[2]);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, argv[2], strlen(argv[2]));
    strncat(command, " \n", 2);

    /* libera memoria */
    free(image_pointer);
    free(image_temp);

    break;

//Filter  Binary
case 3:
    eqMax = 128;

    printf("Filter Binary...\n");

    //initialize
    width_out=width;
    height_out=height;
    bit_depth_out=bit_depth;

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

    for(i=0; i< width_out* height_out; i++)
        image_pointer_out[i] = filterTransform(image_pointer[i]);

    write_pgm_file(argv[2]);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, argv[2], strlen(argv[2]));
    strncat(command, " \n", 2);
    break;

//Filter Log
case 4:
    printf("Filter Log...\n");
    if(argc !=4) abort();

    //initialize
    width_out=width;
    height_out=height;
    bit_depth_out=bit_depth;

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);
    for(i=0; i< width_out* height_out; i++){
          image_pointer_out[i] = filterLog(image_pointer[i]);
    }
    write_pgm_file(argv[2]);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, argv[2], strlen(argv[2]));
    strncat(command, " \n", 2);
    break;

//Filter Exp
case 5:
    printf("Filter Exp...\n");
    if(argc < 3) abort();

    //initialize
    width_out=width;
    height_out=height;
    bit_depth_out=bit_depth;

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);
    for(i=0; i< width_out* height_out; i++){
         image_pointer_out[i] = filterExp(image_pointer[i]);
    }
    write_pgm_file(argv[2]);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, argv[2], strlen(argv[2]));
    strncat(command, " \n", 2);
    break;

case 6:
case 7:
case 8:
    if(argc < 4)abort();
    int type = atoi(argv[3]);
    int dim = atoi(argv[4]);
    if((dim !=3)||(dim !=5)||(dim !=7)){
        if(type !=6){
            printf("ERR | Param incorrect dim: 3, 5, o 7\n\n");
            printf("%i\n",dim);
            //abort();
        }
    }
    char *modeGrad;
    if(dim == 0) modeGrad = "dx+dy";
    if(dim == 1) modeGrad = "dx";
    if(dim == 2) modeGrad = "dy";
    if(type== 6) printf("Apply Filter Gradient: %s \n", modeGrad);
    if(type== 7) printf("Apply Filter Laplacian dim: %dx%d ...\n", dim, dim);
    if(type== 8) printf("Apply Filter Gaussian  dim: %dx%d ...\n", dim, dim);

    if(dim == dim3){
        height_out = height-2;
        width_out = width-2;
    }
    if(dim == dim5){
        height_out = height-4;
        width_out = width-4;
    }
    if(dim == dim7){
        height_out = height-6;
        width_out = width-6;
    }
    bit_depth_out=bit_depth;

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

    if(type == 6)filterGradient(dim);
    if(type == 7)filterLaplacian(dim);
    if(type == 8)filterGaussian(dim);

    write_pgm_file(argv[2]);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, argv[2], strlen(argv[2]));
    strncat(command, " \n", 2);
    break;

/*   Operation change Size proportional argv[4]*/
case 9:
case 10:
case 11:
    if(argc < 4)abort();
    printf("Change Size ...\n");
    int mode = atoi(argv[3]);
    if(mode == 9) printf("Mode: Vecino proximo \n");
    if(mode ==10) printf("Mode: Bilineal \n");
    if(mode ==11) printf("Mode: Bicubic \n");

    float valX, valY;
    bit_depth_out=bit_depth;

    //initialize
    height_out = height*atof(argv[4]);
    width_out = width*atof(argv[4]);

    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

    float ival, jval;
    i=0;
    for(ival=1; ival<=height_out; ival++){
        for(jval=1; jval<= width_out; jval++){
            valX = ival/height_out;
            valY = jval/width_out;
            if(mode ==9) image_pointer_out[i] = inter_Neighbor(valY, valX);
            if(mode ==10)image_pointer_out[i] = inter_lin(valY, valX);
            //if(mode ==11);
            i++;
        }
    }
    write_pgm_file(argv[2]);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, argv[2], strlen(argv[2]));
    strncat(command, " \n", 2);
    break;

//Add Noise Gauss
case 12:
    strcpy(file, "noiseG_");
    strncat(file, argv[1], strlen(argv[1])-4);
    printf("\ngenerate  %s:\n", file);
    printf("%s \n",file);

    height_out = height;
    width_out = width;
    bit_depth_out=bit_depth;
    image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

    double um = atoi(argv[4]);
    double desv =atof(argv[5]);
    printf("um:%.2f, %desv:.2f \n", um, desv);

    double noise;
    for(i=0; i< width_out* height_out; i++){
        noise = generatoNoisedoubleRandn (um, desv);
        if(noise > 0)image_pointer_out[i] = image_pointer[i]+ noise;
        else image_pointer_out[i] = image_pointer[i];
    }
    write_pgm_file(file);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, file, strlen(file));
    strncat(command, " \n", 2);

    break;

case 13:

  //  WriteText();
    convolucionVector();
    break;

//This operation Graphic de histogram
default:
    strcpy(file, "histogram_");
    strncat(file, argv[1], strlen(argv[1])-4);
    printf("\ngraph  %s: \n\n", file);

    //define tam image histogram 256x 256
    width_hist = sizeLevel;
    height_hist = sizeLevel;
    bit_depth_hist=bit_depth;
    image_histogram=(int*)malloc(sizeof(int)*width_hist*height_hist);
    sumHist=(int*)malloc(sizeof(int)*sizeLevel);

    getHistogram(image_pointer);
    //normalizeHistogram();
    //accumHostogram();

    write_pgm_hist(file);

    //open file with imageViewer
    strcpy(command, "eog ");
    strncat(command, file, strlen(file));
    strncat(command, " \n", 2);
    printf(command);
    //system("eog /home/ferminm/build-ProjectImagePro-Desktop-Debug/histogramB\n");
    break;
}

//process_pgm_file();
system(command);
//write_pgm_file(argv[2]);
printf("End program ...\n");

return 0;
}


/* Funcion de Lectura */
void read_pgm_file(char *filename){
int i,c,c0;

/* Abre el archivo de lectura */
FILE *fp = fopen(filename, "r");


/* Lee la  linea inicial y los comentarios */
for(i=0;i<Maxline-1 && (c=fgetc(fp)) != '\n'; i++);
while((c0=fgetc(fp))=='#')
{printf("%c",c0); {for(i=0;i<Maxline-1 && (c=fgetc(fp)) != '\n'; i++) printf("%c",c);} printf("\n");}
ungetc(c0,fp);

/* Lee los parametros de la imagen */
fscanf(fp,"%d",&width);
fscanf(fp,"%d",&height);
fscanf(fp,"%d",&bit_depth);

printf("Parametros de Entrada: w:%d \t h:%d \t b:%d \n",width,height,bit_depth);

image_pointer=(int *)malloc(sizeof(int)*width*height);

for(i=0;i<width*height;i++) fscanf(fp,"%d",&image_pointer[i]);

fclose(fp);
}

/* Funcion de Escritura */
void write_pgm_file(char *filename)
{
int i;

/* Abre el archivo de escritura */
FILE *fp = fopen(filename, "w");

/* Escribe  las dos lineas iniciales */
fprintf(fp,"P2\n");
fprintf(fp,"# Creado por programa xxx\n");


fprintf(fp,"%d\t",width_out);
fprintf(fp,"%d\n",height_out);
fprintf(fp,"%d\n",bit_depth_out);

printf("Parametros de Salida: %d \t %d \t %d \n",width_out,height_out,bit_depth_out);
printf("\n");

for(i=0;i<width_out*height_out;i++) fprintf(fp,"%d\n",image_pointer_out[i]);

fclose(fp);
}

void write_pgm_hist(char *name)
{
int i;

/* Abre el archivo de escritura */
FILE *fp = fopen(name, "w");

/* Escribe  las dos lineas iniciales */
fprintf(fp,"P2\n");
fprintf(fp,"# Histograma\n");

fprintf(fp,"%d\t",width_hist);
fprintf(fp,"%d\n",height_hist);
fprintf(fp,"%d\n",bit_depth_hist);

printf("Parametros de histograma: %d \t %d \t %d \n",width_hist,height_hist,bit_depth_hist);
printf("\n");

for(i=0;i<width_hist*height_hist;i++) fprintf(fp,"%d\n",image_histogram[i]);

fclose(fp);
}


/* Funcion de Procesamiento */
void process_pgm_file()
{
   width_out=width;
   height_out=height;
   bit_depth_out=bit_depth;

   image_pointer_out=(int *)malloc(sizeof(int)*width_out*height_out);

   // Hace algo para calcular image_pointer_out
   //for(i=0;i<width_out*height_out;i++) image_pointer_out[i]=image_pointer[i];
}

void getHistogram(int *inputPointer){
    int i, j, k;
    int valMax=0;
    float sum = 0;
    float sumTot =0;

    //histogram
    for(i=0; i<sizeLevel; i++)histogram[i] = 0;
    for(i=0; i< width*height; i++)histogram[inputPointer[i]] += 1;
    for(i=0; i<sizeLevel; i++){
        printf("%.0f\t", histogram[i]);
        sumTot += histogram[i];
        if(histogram[i]>valMax)valMax = histogram[i]; //valor max
    }
    printf("\n");
    printf("\nMayor is: %d\n", valMax);

    //representacion accumulativa
    for(i=0; i<sizeLevel; i++){
        sum += histogram[i];
        sumHist[i] = sum/ sumTot;
    }
    //for(i=0; i<256;i++)printf("cummHist: %.3f\n", sumHist[i]);

    FILE *fileH = fopen("histogram", "w");
    for(i=0;i<sizeLevel;i++) fprintf(fileH,"%i\t%.1f\n",i,histogram[i]);
    fclose(fileH);

    FILE *fileHA = fopen("histogramAccum", "w");
    for(i=0;i<sizeLevel;i++) fprintf(fileHA,"%i\t%.3f\n",i,sumHist[i]);
    fclose(fileHA);

    //normalización para representation
    for (i=0; i<sizeLevel; i++){
        histogramNormal[i] = histogram[i]/valMax;
        //printf("%03.3f\n", histogramNormal[i]);
    }

    //Grafico de histograma
    k = 0;
    for (j=0; j< sizeLevel; j++){
        i=0;
        while(i< width_hist){
            for(; i<histogramNormal[j]*(sizeLevel-1); i++){
                //printf("*");
                image_histogram[k] =0;
                k++;
            }
            //printf("-");
            i++;
            image_histogram[k] =255;
            k++;
        }
        //printf("\n");
    }
}

/* Funcion  para ecualizacion de imagen - obtiene la distribucion del histograma para luego escalar */
void equalizeImage(){

    float maxGrey=255.0;
    float minGrey=0.0;
    float scale = (maxGrey-minGrey)/(height*width);
    int i, sum;
    unsigned char* accumScale;
    accumScale=(int *)malloc(sizeof(int)*sizeLevel);
    printf("\n");

    i=0;
    sum=0;
    while(i < sizeLevel){
        sum += histogram[i];
        accumScale[i] = (unsigned char)(sum * scale);
        ++i;
    }

    i = 0;
    while(i < (height*width)){
        image_pointer_out[i] = accumScale[image_pointer[i]];
        ++i;
    }
}

/* Funcion de interpolacion Vecino proximo (Neighbor). Las coordenadas deben estar entre 0 y 1 */
float inter_Neighbor(float x,float y)
{
    int ix,iy,i0;
    float a0;

    if(x <0 || x>=1 || y <0 || y>=1) return(0);

    ix= x*width;
    iy= y*height;

    i0=ix+iy*width;
    a0=image_pointer[i0];
    return a0;
}


/* Funcion de interpolacion bi-lineal. Las coordenadas deben estar entre 0 y 1 */
float inter_lin(float x,float y){
    int ix,iy,i0;
    float dx,dy,a0,a1,a2,a3;

    if(x <0 || x>=1 || y <0 || y>=1) return(0);
    //printf("value x: %f, value y: %f\n", x, y);

    ix= x*width;
    iy= y*height;

    dx=x*width-ix;
    dy=y*height-iy;

    i0=ix+iy*width;

    a0=image_pointer[i0];
    a1=image_pointer[i0+1];
    a2=image_pointer[i0+width];
    a3=image_pointer[i0+width+1];

return(a0+(a1-a0)*dx+(a2-a0)*dy+(a0+a3-a1-a2)*dx*dy);
}

/* Funcion de interpolacion bicubica. Las coordenadas deben estar entre 0 y 1 */
float  inter_bicubic(float x, float y){
   float dx,dy,a0,a1,a2,a3, ca,cb;
   int d0,d2,d3;

   d2= x*width;
   d3= y*height;

   dx=x*width-d2;
   dy=y*height-d3;

   d0=d2+d3*width;

   a0=image_pointer[d0];
   a1=image_pointer[d0+1];
   a2=image_pointer[d0+width];
   a3=image_pointer[d0+width+1];

   a0 = dx;
   a1 = -1.0 / 3 * d0 + d2 - 1.0 / 6 * d3;
   a2 = 1.0 / 2 * d0 + 1.0 / 2 * d2;
   a3 = -1.0 / 6 * d0 - 1.0 / 2 * d2 + 1.0 / 6 * d3;
   ca = a0 + a1 * dx + a2 * dx * dx + a3 * dx * dx * dx;


    a0 = dy;
    a1 = -1.0 / 3 * d0 + d2 -1.0 / 6 * d3;
    a2 = 1.0 / 2 * d0 + 1.0 / 2 * d2;
    a3 = -1.0 / 6 * d0 - 1.0 / 2 * d2 + 1.0 / 6 * d3;
    cb = a0 + a1 * dy + a2 *dy*dy + a3 * dy*dy*dy;
    return (ca + cb);
}

/* Filtro Laplaciano */
void filterLaplacian(int dim){

    int filter1[dim3][dim3] = {{  0, -1,  0 },
                               { -1,  4, -1 },
                               {  0, -1,  0 }};

    int filter2[dim5][dim5] = {{ 0, 0,-1, 0, 0},
                               { 0,-1,-2,-1, 0},
                               {-1,-2,17,-2,-1},
                               { 0,-1,-2,-1, 0},
                               { 0, 0,-1, 0, 0}};

    int filter3[dim7][dim7] ={{ 0, 0,-1,-2,-1, 0, 0},
                              { 0,-1,-2,-4,-2,-1, 0},
                              {-1,-2,-4,-8,-4,-2,-1},
                              {-2,-4,-8,100,-8,-4,-2},
                              {-1,-2,-4,-8,-4,-2,-1},
                              { 0,-1,-2,-4,-2,-1, 0},
                              { 0, 0,-1,-2,-1, 0, 0}};

    int j, i;
    int k, l;
    int y=0;
    float pixel_value = 0.0;

    for(k=0; k<(height-(dim-2));k++){
        for(l=0; l<(width-(dim-1));l++){
            for(i=-(dim/2);i<(dim/2 + 1);i++){
                for(j=-(dim/2);j<(dim/2 + 1); j++){
                    if(dim == dim3)pixel_value += filter1[j+(dim/2)][i+(dim/2)]*(float)image_pointer[l+(j+(dim/2))+(k*width)+((i+(dim/2))*width)];
                    if(dim == dim7)pixel_value += filter3[j+(dim/2)][i+(dim/2)]*(float)image_pointer[l+(j+(dim/2))+(k*width)+((i+(dim/2))*width)];
                    if(dim == dim5)pixel_value += filter2[j+(dim/2)][i+(dim/2)]*(float)image_pointer[l+(j+(dim/2))+(k*width)+((i+(dim/2))*width)];
                    //if(k<1)printf("%i\t%i\t%i\t%i\t%i\n",j,i,l,k,l+(j+3)+(k*width)+((i+3)*width));
                }
            }
            if (pixel_value < 0) pixel_value = 0;
            if (pixel_value > 255) pixel_value = 255;
            image_pointer_out[y] = (int)pixel_value;
            pixel_value = 0;
            y++;
       }
    }
}

/* Filtro Gausseano */
void filterGaussian(int dim){

    int filter1[dim3][dim3]= {{1,2,1},
                              {2,4,2},
                              {1,2,1}};

    int filter2[dim5][dim5] ={{1,1,1,1,1},
                              {1,1,1,1,1},
                              {1,1,1,1,1},
                              {1,1,1,1,1},
                              {1,1,1,1,1}};

    int filter3[dim7][dim7] ={{1,1,1,1,1,1,1},
                              {1,1,1,1,1,1,1},
                              {1,1,1,1,1,1,1},
                              {1,1,1,1,1,1,1},
                              {1,1,1,1,1,1,1},
                              {1,1,1,1,1,1,1},
                              {1,1,1,1,1,1,1}};

    int j, i;
    int k, l;
    int y=0;
    float temp = 0.0;
    float pixel_value = 0.0;

    //calculo para normalizacion
    for (k = 0; k <= (dim-1); k++) {
       for (j =0;  j<= (dim-1); j++) {
           if (dim == dim3)temp +=filter1[j][k];
           if (dim == dim5)temp +=filter2[j][k];
           if (dim == dim7)temp +=filter3[j][k];
       }
    }

    for(k=0; k<(height-(dim-2));k++){
        for(l=0; l<(width-(dim-1));l++){
            for(i=-(dim/2);i<(dim/2 + 1);i++){
                for(j=-(dim/2);j<(dim/2 + 1); j++){
                    if(dim == dim3)pixel_value += filter1[j+(dim/2)][i+(dim/2)]/temp * (float)image_pointer[l+(j+(dim/2))+(k*width)+((i+(dim/2))*width)];
                    if(dim == dim5)pixel_value += filter2[j+(dim/2)][i+(dim/2)]/temp * (float)image_pointer[l+(j+(dim/2))+(k*width)+((i+(dim/2))*width)];
                    if(dim == dim7)pixel_value += filter3[j+(dim/2)][i+(dim/2)]/temp * (float)image_pointer[l+(j+(dim/2))+(k*width)+((i+(dim/2))*width)];
                    //if(k<1)printf("%i\t%i\t%i\t%i\t%i\n",j,i,l,k,l+(j+dim/2)+(k*width)+((i+dim/2)*width));
                }
            }
            if (pixel_value < 0) pixel_value = 0;
            if (pixel_value > 255) pixel_value = 255;
            image_pointer_out[y] = (int)pixel_value;
            pixel_value = 0;
            y++;
       }
    }
}

/* Filtro gradiente sqrt(Gx^2+Gy^2) */
void filterGradient(int mode ){

    //Gradiente en x, Gx
    int filter1Gx[dim3][dim3] = {{-1, 0, 1},
                                 {-1, 0, 1},
                                 {-1, 0, 1}};
    //Gradiente en y, Gy
    int filter1Gy[dim3][dim3] = {{ 1, 1, 1},
                                 { 0, 0, 0},
                                 {-1,-1,-1}};
    //G = sqrt(Gx^2 + Gy^2)
    int j, i;
    int k, l;
    int pixel_value;
    int Gx=0;
    int Gy=0;
    int y=0;

    for(l=0; l<(width-1);l++){
        for(k=1; k<(height-1);k++){
            for(i=-1;i<(dim3-1);i++){
                for(j=-1;j<(dim3-1); j++){
                   if((mode ==1)||(mode == 0))Gx += filter1Gx[j+1][i+1]*(float)image_pointer[k+j+(l*width)+((i+1)*height)];
                   if((mode ==2)||(mode == 0))Gy += filter1Gy[j+1][i+1]*(float)image_pointer[k+j+(l*width)+((i+1)*height)];
                }
            }
            pixel_value = sqrt(Gx*Gx + Gy*Gy);
            Gx =0; Gy =0;
            if (pixel_value < 0) pixel_value = 0;
            if (pixel_value > 255) pixel_value = 255;
            image_pointer_out[y] = (int)pixel_value;
            y++;
       }
    }
}

/* Funciones de Transformacion T */
int filterTransform(int valPixel){
    if(valPixel>= eqMax) return 0;
    else return valPixel;
}

float filterLog(int valPixel){
    if(valPixel < 0) valPixel*=(-1);
    return (255/log(1+255))*log(1+valPixel);
}

/* Funtion Exp(c*Pix[ij])*/
float filterExp(int valPixel){
    if(valPixel < 0) valPixel*=(-1);
    return  exp(0.02174*valPixel);
}

/* Generador de ruido gaussiano*/
double generatorNoiseGauss(int m, float d)
{
  double temp1;
  double temp2;
  double result;
  float p;
  // 1/ sqrt(2 * PI * Desv)
  // 0,607 / sqrt(2 * PI * Desv)

  p = d;
  while( p > 0 ){
    temp2 = ( rand() / ( (double)RAND_MAX ) );
    if ( temp2 == 0 ){// temp2 is >= (RAND_MAX / 2)
      p = d;
    }// end if
    else{// temp2 is < (RAND_MAX / 2)
       p = -d;
    }
  }
  temp1 = cos( (2.00 * (double)PI ) * rand()/((double)RAND_MAX));
  result = m*sqrt( -2.00 * log( temp2 )) * temp1;
  return result;
}

double generatoNoisedoubleRandn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1){
      call = !call;
      return (mu + sigma * (double) X2);
  }
  do{
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;
  return (mu + sigma * (double) X1);
}

void WriteTextPoint(){

    printf("TEST\n");
    int i;
    float temp5, temp4;



    FILE *fp3 = fopen("filename", "w");
    for(i=0; i<500; i++){
            temp5 = (cos(i)* rand() / ( (double)RAND_MAX ) );
            temp4 = (sin(i)* rand() / ( (double)RAND_MAX ) );

            printf("%.5f\t %.5f \n",temp4,temp5);
            fprintf(fp3,"%.5f\t %.5f \n",temp4,temp5);

    }
    printf("TEST");

    fclose(fp3);

}

void convolucionVector(){

    //Gradiente en x, Gx
    int filter[dim3][dim3] = {{-1, 0, 1},
                              {-2, 0, 2},
                              {-3, 0, 3}};

    int vector[25] = {1,1,1,1,1,
                      2,2,2,2,2,
                      3,3,3,3,3,
                      4,4,4,4,4,
                      5,5,5,5,5};

    int i,j,k;
    printf("test convolution\n ");

    int l=0;
    k=0;

    for(l=0; l<(dim5-2);l++){
        for(k=1; k<(dim5-1);k++){
            for(i=-1;i<(dim3-1);i++){
                for(j=-1;j<(dim3-1); j++){
                   printf("ij:%i\n", filter[i+1][j+1]);
                   printf("%i\n",k+j+(l*dim5)+((i+1)*dim5));
                }
            }
       }
    }
}
