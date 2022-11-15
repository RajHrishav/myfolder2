#include<stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#define        CHROM   0
#define        POS     1
#define        ID      2
#define        REF     3
#define        ALT     4
#define        QUAL    5
#define        FILTER  6
#define        INFO    7
#define        FORMAT  8
//chr_num vcf_file min_mut max_mut gid_list score_mat num_clust max_iter threshold
//a 20 "sample.vcf" 0 100 "gid.txt" "score_file.dat" 2 100 0.5
float DATA[10000000][50];
int data_row_count;
//float score_mat[4][4];

unsigned long linesToDelete[65000000] ;
int do_validation(char * );
void SORT(int n,char column_names[n][20]);
float map_data(float score_mat[][4],char REFALT_ARR[],char x,char y);
float InsData_calcmutn(char* line_buf,float temp_data[],char gid_list[][20],int gid_list_size,float score_mat[][4],char column_names[][20]);

int main(int argc,char* argv[])
{
		
	

        int chromnbr_ip ;
        char VCF_FILENAME[500];
        float minmutn_ip ;
        float maxmutn_ip ;
		char GID_FILENAME[500];
		char SCORE_FILENAME[500];
		int num_clust;
		int max_iter;
		float threshold;
		int i, j;
		char column_names[3000][20];
		int no_of_vcf_columns;
		char* piece = NULL;
		
		
        char *line_buf = NULL;
        size_t line_buf_size = 0;
        //int line_count = 0;
        ssize_t line_size;
		char tmpCharptr1[20000];
		int is_valid;
		int chromnbr;
		float mutn;
		unsigned long result_count;
		result_count=0;
		
		printf("\n60\n");
		
        chromnbr_ip = atoi(argv[1]);
        strcpy(VCF_FILENAME,argv[2]);
        minmutn_ip = atof(argv[3]);
        maxmutn_ip = atof(argv[4]);
		strcpy(GID_FILENAME,argv[5]);
		strcpy(SCORE_FILENAME,argv[6]);
		num_clust = atoi(argv[7]);
		max_iter = atoi(argv[8]);
		threshold = atof(argv[9]);
		
		FILE *fp;
		char gid_list[2500][20];
		int gid_list_size = 0;
		float temp_data[gid_list_size];//will store corresponding mapped data for the columns mentioned in GID file for 1 row
		float score_mat[4][4];
		printf("\n76\n");
		
		fp = fopen(GID_FILENAME, "rb");
		line_size = getline(&line_buf, &line_buf_size, fp);
		while (line_size >= 0) {
			strcpy( gid_list[gid_list_size++],line_buf );
			line_size = getline(&line_buf, &line_buf_size, fp);
		}
		free(line_buf);
		line_buf = NULL;
		fclose(fp);
		fp = NULL;
			
		//printf("\n90\n%s",SCORE_FILENAME);	
		//reading data from score file
		fp = fopen(SCORE_FILENAME,"rb");
		if(fp) {
			for(i = 0; i <= 3; i++) {
				for(j = 0; j <= 3;j++) {
					if(!fscanf(fp,"%f",&score_mat[i][j])) break;
				}
			}
			fclose(fp);
			fp = NULL;
			free(line_buf);
			line_buf = NULL;
		}
		
	fp = fopen(VCF_FILENAME, "rb");
        if (!fp)        {
                fprintf(stderr, "Error opening file '%s'\n", VCF_FILENAME);
                return 1;
        }
		
        line_size = getline(&line_buf, &line_buf_size, fp);
		//skipping lines starting with ##
        while(*line_buf == '#' && *(line_buf+1) == '#'){
            line_size = getline(&line_buf, &line_buf_size, fp);
        }
		
		//saving column names in an array from the line starting with single #
	if(*line_buf == '#' && *(line_buf+1) != '#'){
			piece = strtok(line_buf,"	");
			no_of_vcf_columns = 0;
			while(piece != NULL) {
				strcpy(column_names[no_of_vcf_columns],piece);			
				piece = strtok(NULL,"	");
				no_of_vcf_columns++;
			}
			
			//SORT(no_of_vcf_columns,column_names);
            line_size = getline(&line_buf, &line_buf_size, fp);
        }
        
        printf("\n128\n");
	//processing actual data
        while (line_size >= 0) {
                strcpy(tmpCharptr1,line_buf);
                //printf("\n%s",tmpCharptr1);
                is_valid = do_validation(tmpCharptr1);
                printf("\n%d",is_valid);
                if(is_valid == 1)	{
                    piece = strtok(tmpCharptr1,"    ");
                    chromnbr = atoi(piece);
                    if(chromnbr == chromnbr_ip){
                        //calculating mutation & inserting gid columns data in a temporary array
						mutn = InsData_calcmutn(line_buf,temp_data,gid_list,gid_list_size,score_mat,column_names);
                        printf("\nmutn%f",mutn);
                        if((mutn >= minmutn_ip) && (mutn <= maxmutn_ip)){
							
                            result_count++;//increasing valid row count
							
							for(j = 0; j <= gid_list_size-1; j++) { // inserting data in original data matrix from temporay array
								DATA[data_row_count][j] = temp_data[j];
							}
							
							data_row_count++;
                        }
                    }
                }
			line_size = getline(&line_buf, &line_buf_size, fp);
        }//1st while i.e while(line_size>0)
        free(line_buf);
        line_buf = NULL;
        fclose(fp);
        
        printf("\n158\n");
        printf("%lu\n",result_count);
        printf("%d %d",data_row_count,gid_list_size);
		
	for(i = 0; i <= data_row_count-1; i++) {
		for(j = 0; j <= gid_list_size-1; j++) {
			//printf("%f\t",DATA[i][j]);
			printf("%d-%d",i,j);
		}
		printf("\n");
	}

		
		printf("hi");
        return 0;
}//main method


int do_validation(char * tmpCharptr)
{
	static unsigned long line_count = 0;
	static unsigned long i= 0;
	int tab_count, j, flg;

	line_count++;
	tab_count = 0;

    if(*tmpCharptr=='#'){
        linesToDelete[i++]=line_count;//storing line numbers that needs to be deleted.
			return 0;
    }
	else{
			//chrom column Integer check
			if(tab_count==CHROM){
					flg=0;
					while(*tmpCharptr!='\t'){
							if(*tmpCharptr<'0' ||  *tmpCharptr>'9'){
									linesToDelete[i++]=line_count;
									flg=1;
									printf("\n201-%c",*tmpCharptr);
									//break;
									return 0;
							}
							tmpCharptr++;
					}
					if(flg==0){
							while(tab_count!=REF){
									if(*tmpCharptr=='\t'){
											tab_count++;
									}
									tmpCharptr++;
							}
					}
			}
			if(tab_count==REF){
					flg=0;
					while(*tmpCharptr!='\t'){
							//REF not belongs to set {A,T,G,C}
							if((*tmpCharptr!='A') && (*tmpCharptr!='T') && (*tmpCharptr!='G') && (*tmpCharptr!='C')){
									linesToDelete[i++]=line_count;
									flg=1;
									printf("\n223");
									//break;
									return 0;
							}
							//more than one character in REF
							else if(*(tmpCharptr+1)!='\t'){
									linesToDelete[i++]=line_count;
									flg=1;
									printf("\n231");
									//break;
									return 0;
							}
							tmpCharptr++;
					}
					if(flg==0){
							while(tab_count!=ALT){
									if(*tmpCharptr=='\t'){
											tab_count++;
									}
									tmpCharptr++;
							}
					}
			}
			if(tab_count==ALT){
					flg=0;
					j=1;
					while(*tmpCharptr!='\t'){
							if((j%2==1) && ((*tmpCharptr!='A') && (*tmpCharptr!='T') && (*tmpCharptr!='G') && (*tmpCharptr!='C'))){
									linesToDelete[i++]=line_count;
                                    flg=1;
									//break;
									printf("\n254");
									return 0;
							}
							else if((j%2==0) && (*tmpCharptr!=',')){
									linesToDelete[i++]=line_count;
									flg=1;
									//break;
									printf("\n261");
									return 0;
							}
							else if((j%2==0) && (*(tmpCharptr+1)=='\t') && (*tmpCharptr==',')){
											linesToDelete[i++]=line_count;
											flg=1;
											//break;
											printf("\n268");
											return 0;
											}
							j++;
							tmpCharptr++;
					}
					if(flg==0){
							while(tab_count!=FILTER){
									if(*tmpCharptr=='\t'){
											tab_count++;
									}
									tmpCharptr++;
							}

					}
			}

			if(tab_count==FILTER){
					if(*tmpCharptr != 'P'){
						linesToDelete[i++]=line_count;
						printf("\n288");
						return 0;
					}
			}

	}//1st else
	return 1;
}

float map_data(float score_mat[][4],char REFALT_ARR[],char x,char y)
{
	int indx_i,indx_j;
	x = REFALT_ARR[x - '0'];
	y = REFALT_ARR[y - '0'];
	switch(x)
	{
		case 'A' :
			indx_i = 0;
			break;
		case 'C' :
			indx_i = 1;
			break;
		case 'G' :
			indx_i = 2;
			break;
		case 'T' :
			indx_i = 3;
			break;
	} 
	
	switch(y)
	{
		case 'A' :
			indx_j = 0;
			break;
		case 'C' :
			indx_j = 1;
			break;
		case 'G' :
			indx_j = 2;
			break;
		case 'T' :
			indx_j = 3;
			break;
	}
	
	return score_mat[indx_i][indx_j]; // map_matrix is 4x4 matrix .Will be given as input.
}

float InsData_calcmutn(char* line_buf,float temp_data[],char gid_list[][20],int gid_list_size,float score_mat[][4],char column_names[][20])
{
	int j,gid_indx,r;
	for(j=0; j <= gid_list_size-1; j++) { //reset before insertion
		temp_data[j] = 0;
	}
	
	int tab_count, mutn00_count,totalmutn_count;
	char REFALT_ARR[55];
	char tmpCharptr[20000];
	
	strcpy(tmpCharptr,line_buf);
	mutn00_count = 0;
	totalmutn_count = 0;
	tab_count = 0;
	gid_indx = 0;
	r = 0;
	char * piece = strtok(tmpCharptr,"	");
	while(piece != NULL) {
		
		if(tab_count > FORMAT){// format means 8
			if(strlen(piece) >= 3 && *(piece + 1) == '|'){
				if(*piece == '0' && *(piece + 2) == '0'){
					mutn00_count++;
				}
				totalmutn_count++;
			}
		}
		
		if(tab_count == REF){REFALT_ARR[r++] = piece[0];}
		if(tab_count == ALT){while(*piece != '\0') { REFALT_ARR[r++] = *piece; piece++;}};//insert alt in refalt_arr char array
		
		if(strcmp(column_names[tab_count],gid_list[gid_indx]) == 0){
			if(strlen(piece) >= 3 && *(piece + 1) == '|'){
				temp_data[gid_indx] = map_data(score_mat,REFALT_ARR,piece[0],piece[2]);
				gid_indx++;
			}
		}

		piece = strtok(NULL,"	");
		tab_count++;
	}
	
	return ((float)(totalmutn_count - mutn00_count)/totalmutn_count) * 100;
}

void SORT(int n,char array[n][20])
{
  char temp[20];
  for(int i=0; i<=n-1; i++){
    for(int j=0; j< 20; j++){
      if(strcmp(array[j], array[j+1]) > 0){
        //swap array[j] and array[j+1]
        strcpy(temp, array[j]);
        strcpy(array[j], array[j+1]);
        strcpy(array[j+1], temp);
      }
    }
  }
}
