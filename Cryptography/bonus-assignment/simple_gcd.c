#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <math.h>

int count_lines_in_file(FILE* fptr);
void read_numbers(mpz_t arr[], int n, FILE* fptr);
int build_product_tree(mpz_t arr[], int tier_count, int n,  int tier_lengths[]);
int fill_product_tree_tier(mpz_t arr[], int tier_count, int tier, int start, int end, int tier_lengths[]);
void transform_to_remainder_tree(mpz_t arr[], int tier, int start, int prev_start, int tier_lengths[]);
void find_factors(mpz_t remainder_tree[], mpz_t arr[], int n);
void factor_outliers(mpz_t arr[], int outliers[], int outlier_count, int factorized[], int factorized_count);


int main(int argc, char const *argv[]) {
  time_t begin = time(NULL);  

  // read info from file
  FILE* fptr = fopen("moduli.txt", "r");
  int n = count_lines_in_file(fptr);
  mpz_t arr[n]; //load numbers from file in array
  read_numbers(arr, n, fptr);
  fclose(fptr);

  //build product tree
  mpz_t product_tree[2*n+20];
  for (int i = 0; i < n; i++) {
    mpz_init_set(product_tree[i], arr[i]); //copy elements from arr to product tree
  }
  int tier_count = (int)ceil(log2((double)n)) + 1;
  int tier_lengths[tier_count]; //to remember lengths of each tier
  int end_index = build_product_tree(product_tree, tier_count, n, tier_lengths);
  time_t end_0 = time(NULL);
  printf("--------------------------------\n");
  printf("Building product tree took %ld seconds\n", (end_0 - begin));
  printf("--------------------------------\n");
  
  //transform product tree to remainder tree
  transform_to_remainder_tree(product_tree, tier_count-2, end_index-1, end_index, tier_lengths);
  time_t end_1 = time(NULL);
  printf("--------------------------------\n");
  printf("Transforming product tree to remainder tree took %ld seconds\n", (end_1 - end_0));
  printf("--------------------------------\n");

  //find factors and write them to file
  find_factors(product_tree, arr, n); //product tree here is remainder tree
  time_t end_2 = time(NULL);
  printf("Finding factors took %ld seconds\n", (end_2 - end_1));
   printf("--------------------------------\n");

  //free memory
  for (int i = 0; i < n; i++) {
    mpz_clear(arr[i]);
  }
  for (int i = 0; i < 2*n+20; i++) {
    mpz_clear(product_tree[i]);
  }

  time_t end = time(NULL);
  printf("Whole run took %ld seconds\n", (end - begin));
  printf("--------------------------------\n");
  return 0;
}

int count_lines_in_file(FILE* fptr){
  char buffer[1024];
  // count lines in file
  int n = 0;
  while (fgets(buffer, 1024, fptr) != NULL) {
    n++;
  }
  // rewind file
  rewind(fptr);
  return n;
}

void read_numbers(mpz_t arr[], int n, FILE* fptr){
  // read n lines as array of gmp integer objects from file
  char buffer[1024];
  for (int i = 0; i < n; i++) {
    fscanf(fptr, "%1023s", buffer);
    mpz_init_set_str(arr[i], buffer, 16);
  }
}

int build_product_tree(mpz_t arr[], int tier_count, int n,  int tier_lengths[]){
  //initialize mpz objects from n to end of arr
  for (int i = n; i < 2*n + 20; i++) {
    mpz_init(arr[i]);
  }
  // build product tree
  tier_lengths[0] = n;
  int end_index = fill_product_tree_tier(arr, tier_count, 2, 0, n, tier_lengths);

  return end_index;
}

int fill_product_tree_tier(mpz_t arr[], int tier_count, int tier, int start, int end, int tier_lengths[]){
  int new_index = end;
  for (int i = start; i < end; i++){
    if (i + 1 >= end){
      mpz_set(arr[new_index], arr[i]);
    }else{
      mpz_mul(arr[new_index], arr[i], arr[i+1]);
      i++;
    }
    new_index++;
  }
  int end_index = new_index-1;
  tier_lengths[tier-1] = new_index - end;

  printf("Filled product tree tier %d, Indices: %d -- %d\n", tier, end, end_index);
  if (tier < tier_count){
    end_index = fill_product_tree_tier(arr, tier_count, tier+1, end, new_index, tier_lengths);
  }

  return end_index;
}


void transform_to_remainder_tree(mpz_t arr[], int tier, int start, int prev_start, int tier_lengths[]){
  double prev_tier_index = (double)prev_start;
  if (tier_lengths[tier] % 2 == 1){
    prev_tier_index -= 0.5;
  }

  int end = start - tier_lengths[tier];
  for (int i = start; i > end; i--){
    //calculate square of arr[i]
    mpz_t square;
    mpz_init(square);
    mpz_mul(square, arr[i], arr[i]);
    // set arr[i] to arr[prev_tier_index] modulo square
    int idx = (int)ceil(prev_tier_index);
    mpz_mod(arr[i], arr[idx], square);
    prev_tier_index -= 0.5;
    if (prev_tier_index < prev_start - tier_lengths[tier+1]){
      break;
    }
  }

  printf("Filled remainder tree tier %d, Indices: %d -- %d\n", tier+1, start, end+1);
  if (tier > 0){
    transform_to_remainder_tree(arr, tier-1, end, start, tier_lengths);
  }
}


void find_factors(mpz_t remainder_tree[], mpz_t arr[], int n){
  int outliers[n];
  int factorized[n];
  int factorized_count = 0;
  int outlier_count = 0;

  for (int i = 0; i < n ; i++) {
    mpz_t square;
    mpz_init(square);
    mpz_mul(square, arr[i], arr[i]);

    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, square, remainder_tree[i]);

    if (mpz_cmp(gcd, square) == 0){
      outliers[outlier_count] = i;
      mpz_clear(square);
      mpz_clear(gcd);
      continue;
    }    

    // if gcd is more than arr[i] then arr[i] can be factored
    if (mpz_cmp(gcd, arr[i]) > 0){
      printf("FOUND FACTORIZATION OF %d'th number\n", i);
      factorized[factorized_count] = i;
      factorized_count++;

      mpz_t x; // x is gcd/arr[i]
      mpz_init(x);
      mpz_divexact(x, gcd, arr[i]);

      //compute gcd of x and arr[i]
      mpz_t gcd_1;
      mpz_init(gcd_1);
      mpz_gcd(gcd_1, x, arr[i]);

      mpz_t factor;
      mpz_init(factor);
      mpz_divexact(factor, arr[i], gcd_1);

      //append i, factor and gcd to file as new lines
      FILE *fptr;
      fptr = fopen("factors.txt", "a");
      fprintf(fptr, "\n%d\n", i);
      mpz_out_str(fptr, 16, factor);
      fprintf(fptr, "\n");
      mpz_out_str(fptr, 16, gcd_1);
      fclose(fptr);

      // free some memory
      mpz_clear(factor);
      mpz_clear(gcd_1);
      mpz_clear(x);
    }

    // free some memory
    mpz_clear(square);
    mpz_clear(gcd);
  }
  factor_outliers(arr, outliers, outlier_count, factorized, factorized_count);
}

void factor_outliers(mpz_t arr[], int outliers[], int outlier_count, int factorized[], int factorized_count){
   for (int i = 0; i < outlier_count; i++){
    for (int j = 0; j < factorized_count; j++){
        mpz_t gcd;
        mpz_init(gcd); //initialize gmp integer object to 0
        mpz_gcd(gcd, arr[outliers[i]], arr[factorized[j]]); // calculate gcd of i and j
        if (mpz_cmp_ui(gcd, 1) == 0) { // if gcd is 1 arr[i] = gcd * factor1, arr[j] = gcd * factor2
        
          mpz_t factor1;
          mpz_init(factor1);
          mpz_divexact(factor1, arr[outliers[i]], gcd); // factor1 = arr[i] / gcd

          mpz_t factor2;
          mpz_init(factor2);
          mpz_divexact(factor2, arr[factorized[j]], gcd); // factor2 = arr[j] / gcd

          FILE *fptr;
          fptr = fopen("factors.txt", "a"); 
          
          fprintf(fptr, "\ni: %d, j: %d\n", i, j); // write i and j in myfile
          
          mpz_out_str(fptr, 16, factor1); // write factor1 to myfile
          fprintf(fptr, "\n");
          
          mpz_out_str(fptr, 16, gcd); // write gcd to myfile
          fprintf(fptr, "\n");
          
          mpz_out_str(fptr, 16, factor2); // write factor2 to myfile

          // free some memory
          fclose(fptr);
          mpz_clear(factor1);
          mpz_clear(factor2);
        }
        printf("DONE: %d %d\n", i, j); // print which numbers are processing rn
        mpz_clear(gcd);  // free some memory
    }
  }
}