#include "Matmod2.h"

int n;
int k;
int w;
Matmod2 H;
Matmod2 s;

bool init(){
    string filename("input.txt");
    vector<string> lines;
    string line;

    ifstream input_file(filename);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << filename << "'" << endl;
        return false;
    }

    while (getline(input_file, line)){
        lines.push_back(line);
    }

    input_file.close();

    n = stoi(lines[1]); 
    w = stoi(lines[5]); 
    k = n - lines[7].size();

    char c;
    int i = 7;
    vector<string> H_str;
    vector<string> s_str;
    while(true){
        H_str.push_back(lines[i]);
        i++;
        c = lines[i][0];
        if(c == '#') break;
    }

    s_str.push_back(lines[i+1]);
    s = Matmod2(s_str);

    auto extra = Matmod2(H_str);
    H = Matmod2(n-k);
    H.concat(extra);

    return true;    
}

//gauss with only first p columns, and check their rank.
bool gauss(Matmod2& A, Matmod2& b, int p){
    bool is_zero = true;
    for(int t = 0; t < p; t++){
        for(int i = 0; i < n - k; i++){
            if(A.mat[i][t]){
                swap(A.mat[i], A.mat[t]);
                swap(b.mat[i], b.mat[t]);
                for(int j = i + 1; j < n; j++){
                    if(A.mat[j][t]){
                        for(int x = 0; x < n; x++){
                            A.mat[j][x] = (A.mat[j][x] != A.mat[t][x]);
                        }
                        b.mat[j][0] = (b.mat[j][0] != b.mat[t][0]);
                    }
                }
                is_zero = false;
                break;
            }
        }
        if(is_zero){
            return false;
        }
    }

    return true;
}
 
int main()
{
    /* input matrix */
    init();
 
    H.print();
    s.print();

    /*
    gauss(H, s, 8);
    H.print();
    s.print();
    */

 
    return 0;
}