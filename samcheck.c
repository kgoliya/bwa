#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "kseq.h"

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <poll.h>
#include <sys/time.h>


#define	MEM_16G		(1ULL << 34)
#define QUEUE_SIZE 100
#define INIT_OVERFLOW_TABLE_SIZE 1000000
KSEQ_DECLARE(gzFile)


int sam_verbose = 1;

void samcheck_usage(){
    fprintf(stderr,"./bwa samcheck -I1 <first sam> -I2 <second sam> -O <output file> -v <verbose level>\n");
    fprintf(stderr,"Verbose level : 1 -> Only prints how many lines are different\n");
    fprintf(stderr,"Verbose level : 2 -> Prints which lines are different and where\n");
    fprintf(stderr,"-O and -v flags are optional. Default output is stdout and default verbose level is 1\n");
    fprintf(stderr,"This code is only comparing read name and flags. It will ignore header files\n");
}

// Returns -1 when sam record in line1 and line2 are not the same
// Returns 0 when sam record in line1 and line2 are the same
int check_record(char * line1, size_t line1_size,char * line2, size_t line2_size ){
    int i = 0;
    int field_num = 1;
    int lsize = line1_size < line2_size ? line1_size : line2_size;

    for(i = 0;i<lsize;i++){
        if(line1[i] == '\t'){
            field_num++;
            if(field_num >= 3){
                break;
            }
        }
        else{
            //check characters
            if(line1[i] != line2[i]){
                return -1;    
            }
        }

    }
    
    return 0;
}

void check_sam_files(char * sam1_file_name, char * sam2_file_name, char * out_file_name, int output_file){
    FILE * sam1;
    FILE * sam2;
    FILE * sout;

    sam1 = fopen(sam1_file_name,"r");
    sam2 = fopen(sam2_file_name,"r");
    if(output_file == 1){
        sout = fopen(out_file_name,"w");
    }
    else{
        sout = stdout;
    }
    int ret;
    char * line1 = NULL;
    size_t line1_size = 0;

    char * line2 = NULL;
    size_t line2_size = 0;

    int ch = 0;
    int64_t num_diff_lines = 0;
    int64_t line_num = 0;

    while(!feof(sam1) && !feof(sam2)){
        ch = 0;
        line_num++;
        ret = getline(&line1,&line1_size,sam1);
        ret = getline(&line2,&line2_size,sam2);
        
        ch = check_record(line1,line1_size, line2,line2_size);
        if(ch == -1){
            num_diff_lines++;
            if(sam_verbose >= 2){
                fprintf(sout,"[1,%ld] %s",line_num,line1);
                fprintf(sout,"[2,%ld] %s",line_num,line2);
                fprintf(sout,"--------------\n");
            }
        }
    }

    fprintf(stderr,"Number of lines that differ : %ld\n",num_diff_lines);
    if(sam_verbose >= 2){
        if(feof(sam1) && !feof(sam2)){
            fprintf(sout,"%s is smaller",sam1_file_name);
        }
        else if(!feof(sam1) && feof(sam2)){
            fprintf(sout,"%s is smaller",sam2_file_name);
        }
    }


    fclose(sam1);
    fclose(sam2);
    if(output_file == 1){
        fclose(sout);
    }
    return;
}

int main_samcheck(int argc, char *argv[]){
    
    int rargc = 4;
    int i = 0;
    int output_file = 0; 
    int file1 = 0, file2 = 0;
    char * sam1_file_name;
    char * sam2_file_name;
    char * out_file_name;


    if(argc < rargc){
        samcheck_usage();
        return 0;
    }
   
    for(i=0;i<argc;i+=2){
        if(strcmp(argv[i],"-I1") == 0) {
            if(argc > i+1){
                sam1_file_name = argv[i+1];
                file1 = 1;
            }
            else{
                samcheck_usage();
                return -1;
            }
            continue;
        }
        if(strcmp(argv[i],"-I2") == 0) {
            if(argc > i+1){
                sam2_file_name = argv[i+1];
                file2 = 1;
            }
            else{
                samcheck_usage();
                return -1;
            }
            continue;
        }
        if(strcmp(argv[i],"-O") == 0) {
            if(argc > i+1){
                out_file_name = argv[i+1];
                output_file = 1;
            }
            else{
                samcheck_usage();
                return -1;
            }
            continue;
        }
        if(strcmp(argv[i],"-v") == 0){
            if(argc > i+1){
                sam_verbose = (int)strtol(argv[i+1],NULL,10);
                if(sam_verbose < 1){
                    sam_verbose = 1;
                }
                else if(sam_verbose > 2){
                    sam_verbose = 2;
                }
            }
            else{
                samcheck_usage();
                return -1;
            }
            continue;
        }
    }


    check_sam_files(sam1_file_name, sam2_file_name, out_file_name, output_file);


    //fclose(temp_unsorted_file);
    return 0;

}

