#include<time.h>
#include<math.h>

#include "orthros.c"




double checkdl(int argc, char *argv[]) {
	// Define the input/output masks and other parameters for linear analysis
	//#######################################################################
	if (argc < 2) {
		printf("No input provided. Using default values and running the program.\n");
	} else if (argc != 9) {
		printf("Usage: %s <num_of_rounds> <rndOffset> <inputdiff0> <inputdiff1> <outputmask> <DEG> <N> <mode>\n", argv[0]);
		return 1;
	}

	int num_of_rounds = (argc > 1) ? atoi(argv[1]) : 5;
	int rndOffset = (argc > 2) ? atoi(argv[2]) : 4;
	const char *inputdiff0 = (argc > 3) ? argv[3] : "00040000000000000000000000000000";
	const char *inputdiff1 = (argc > 4) ? argv[4] : "00040000000000000000000000000000";
	const char *outputmask = (argc > 5) ? argv[5] : "88800000000000000000000000000000";
	int DEG = (argc > 6) ? atoi(argv[6]) : 20;
	int N = (argc > 7) ? atoi(argv[7]) : 1;
	int mode = (argc > 8) ? atoi(argv[8]) : 1; // 0 for left, 1 for right, 2 for prf
	//#######################################################################

	// Print configuration
	printf("=============================================================\n");
	printf("Configuration:\n");
	printf("  Rounds:         %d\n", num_of_rounds);
	printf("  Offset:         %d\n", rndOffset);
	printf("  Input Diff 0:   %s\n", inputdiff0);
	printf("  Input Diff 1:   %s\n", inputdiff1);
	printf("  Output Mask:    %s\n", outputmask);
	printf("  Degree (DEG):   %d (N1 = 2^%d = %llu samples per experiment)\n", DEG, DEG, 1ULL << DEG);
	printf("  Experiments:    %d\n", N);
	const char *mode_str = (mode == 0) ? "Left Branch" : (mode == 1) ? "Right Branch" : "PRF (both branches)";
	printf("  Mode:           %d (%s)\n", mode, mode_str);
	printf("=============================================================\n");
	printf("Running correlation estimation...\n\n");

	uint64_t N1 = 1ULL << DEG; // Number of queries:  N1 = 2^(DEG)
	int64_t CORR;
	int64_t sum = 0;

	// clock_t clock_timer;
	// clock_timer = clock();

	unsigned char inpDiff0[32],inpDiff1[32], outMask[32];
	hex_string_to_bytes((char*)inputdiff0, inpDiff0, sizeof(inpDiff0));
	hex_string_to_bytes((char*)inputdiff1, inpDiff1, sizeof(inpDiff1));
	hex_string_to_bytes((char*)outputmask, outMask, sizeof(outMask));

	if((mode==modeLeftBr) && (isZero(inpDiff0,sizeof(inpDiff0)))){
		printf("Error: Set inputdiff0");
		return 1;
	} else if((mode==modeRightBr) && (isZero(inpDiff1, sizeof(inpDiff1)))){
		printf("Error: Set inputdiff1");
		return 1;
	}

	for(int num_of_experiment=0;num_of_experiment<N;num_of_experiment++){
		init_prng(num_of_experiment);
		uint64_t counter0 = 0ULL;
		uint64_t counter1 = 0ULL; 
		unsigned char plaintext1[32], plaintext2[32], plaintext3[32];
		unsigned char key[32];
		unsigned char ciphertext1[32], ciphertext2[32];
		
		for(int i=0;i<32;i++){
			key[i]=(unsigned char)(rand() & 0xf);
		}

		for (uint64_t loopcnt = 0; loopcnt < N1; loopcnt++){
			for(int i=0;i<32;i++){
				plaintext1[i]=(unsigned char)(rand() & 0xf);
				plaintext2[i]=plaintext1[i]^inpDiff0[i];
				plaintext3[i]=plaintext1[i]^inpDiff1[i];
			}
			ORTHROS(plaintext1, ciphertext1, key, rndOffset, num_of_rounds, mode, plaintext1);
			ORTHROS(plaintext2, ciphertext2, key, rndOffset, num_of_rounds, mode, plaintext3);
			if(dot_prod(ciphertext1,outMask,sizeof(outMask)) == dot_prod(ciphertext2,outMask,sizeof(outMask))){
				counter0+=1;
			}
			else{
				counter1+=1;
			}   
							
		}
		CORR = (int64_t)(counter0 - counter1);
		sum += CORR;
	}	
	double corr = 0;
	corr = ((double)sum / N) / (1ULL << DEG); // Compute correlation directly without logarithm
	return corr;
}

int main(int argc, char *argv[]) {
	double result = checkdl(argc, argv);
	double log_result = (result == 0) ? 0 : log2(fabs(result));
	int sign = (result > 0) ? 1 : -1;
	printf("=============================================================\n");
	printf("Results:\n");
	printf("  Correlation:                 %lf\n", result);
	printf("  Log2(|Correlation|):         %lf\n", log_result);
	printf("  Sign:                        %d\n", sign);
	printf("=============================================================\n");
	return 0;
}
