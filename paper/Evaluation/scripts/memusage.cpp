#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
using namespace std;
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

string getAbsolutePath(const char* file, char* path, int mode) {
	char* tok;

	/*parse the PATH enviroment*/
	tok = strtok((char*) path, ":");
	while (tok) {
		/*form a new absolute path*/
		string filePath = tok;
		filePath += "/";
		filePath += file;

		/*test the file attribute*/
		if (access(filePath.c_str(), mode) == 0) {
			return filePath;
		}
		/*get the next*/
		tok = strtok(NULL, ":");
	}
	return string("");
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cerr << endl << "MEMUSAGE (by Liu Yongchao, liuy@uni-mainz.de) measures the peak virtual and resident memory of a program" << endl
				 << "--------------------------------------------------------------------------------------------------------" << endl
		 		 <<"Please specify parameters: program [program argument list]" << endl;
		return -1;
	}
	/*get the command name*/
	string command = argv[1];

	/*enviroment variables*/
	string path = getenv("PATH");
	/*check the existence of the command file*/
	if (access(command.c_str(), X_OK) != 0) {
		string newPath = getAbsolutePath((char*) command.c_str(),
				(char*) path.c_str(), X_OK);
		if (newPath.length() > 0) {
			command = newPath;
		} else {
			cerr << "Failed to find the specified command: " << command << endl;
			return -1;
		}
	}
	vector<char*>arglist;
	for(int i = 1; i < argc; ++i){
		arglist.push_back(strdup(argv[i]));
	}
	arglist.push_back(NULL);

	/*create a thread*/
	pid_t pid = fork();
	if (pid == 0) { /*child process*/
		if (execvp(command.c_str(), arglist.data()) == -1) {
			cerr << "Failed to launch the command: " << command << endl;
		}
		return 0;
	} else { /*parent process*/
		if (pid < 0) {
			cerr
					<< "Process creating for the external command failed and will use internal command instead"
					<< endl;
		} else if (pid > 0) {
			int status;
			long peakResidentMem = 0, peakVirtualMem = 0;
			char buffer[1024];
			sprintf(buffer, "ps -p %d -o vsize=,rss=", pid);
			string command = buffer; /*get the command line*/

			/*get the output file name*/
			sprintf(buffer, ".memusage%d", getpid());
			string fileName = buffer;

			/*form the command line*/
			command += " > ";
			command += fileName;
			cerr << "command: " << command << endl;

			/*start the main loop*/
			while (1) {
				/*wait for a time*/
				usleep(100);
				
				/*check the exit of child process*/
				waitpid(pid, &status, WNOHANG);

				/*get the resident memory size*/
				system(command.c_str());

				/*read the first line from the output file*/
				FILE* fp = fopen(fileName.c_str(), "rb");
				if (!fp) {
					cerr << "Failed to open file " << fileName << endl;
					exit(0);
				}
				/*read the file line*/
				if (fgets(buffer, 1023, fp) == NULL) {
					cerr << "The child process exited" << endl;
					break;
				}
				/*close the file*/
				fclose(fp);

				/*parse the line*/
				for (int i = strlen(buffer) - 1;
						i >= 0 && (i == '\n' || i == '\r'); --i) {
					buffer[i] = '\0';
				}
				if (strlen(buffer) == 0) {
					cerr << "The child process exited with buffer length " << strlen(buffer) << endl;
					break;
				}
				/*parse the buffer*/
				static const char delimiter[] = "\t ";
				char* tok = strtok(buffer, delimiter);
				/*get the virtual memory size*/
				if(!tok){
					cerr << "It seems that there is an error in the file: " << fileName << endl;
					break;
				}
				long value = atol(tok);
				peakVirtualMem = max(peakVirtualMem, value);
				/*get the resident memory*/
				tok = strtok(NULL, delimiter);
				if(!tok){
					cerr << "It seems that there is an error in the file: " << fileName << endl;
					break;
				}
				value = atol(tok);
				peakResidentMem = max(peakResidentMem, value);
				
			}

			/*output the peak resident memory*/
			cerr << "Peak resident memory is " << peakResidentMem << " KB; " << (double)peakResidentMem / 1024 << " MB; " << (double)peakResidentMem / 1024 / 1024 << " GB" << endl;

			cerr << "Peak virtual memory is " << peakVirtualMem << " KB; " << (double)peakVirtualMem / 1024 << " MB; " << (double)peakVirtualMem / 1024 / 1024 << " GB" << endl;

			/*remove the file*/
			remove(fileName.c_str());

			return 0;
		}
	}
}

