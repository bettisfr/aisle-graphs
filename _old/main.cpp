// g++ -O3 main.cpp -o oprsc
//
// oprsc.exe type [input] budget rows columns
//
// Example: oprsc.exe 1 32 130 10 20 <--- v1
// Example: oprsc.exe 2 0 130 10 20 <--- v2

#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

void run_oprsc_v1(string input, int B, int m, int n);
void run_oprsc_v2(int B, int m, int n);

int main(int argc, char *argv[]) {
    
    int type = stoi(argv[1]);
    string input = argv[2];
    int B = stoi(argv[3]);
    int m = stoi(argv[4]);
    int n = stoi(argv[5]);
    
    if (type == 1) {
        run_oprsc_v1(input, B, m, n);
    } else if (type == 2) {
        run_oprsc_v2(B, m, n);
    }
}

void run_oprsc_v2(int B, int m, int n) {
    float **r = new float*[m];
    for (int i = 0; i < m; i++) {
        r[i] = new float[n];
    }
    
    string filename = "input/tmp-G.csv";
    ifstream data(filename);
    string line;
    int i = 0;
    while(getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        int j = 0;
        while(getline(lineStream, cell, ',')) {
            r[i][j++] = stoi(cell);
        }
        i++;
    }
    
    float **T = new float*[m];
    for (int i = 0; i < m; i++) {
        T[i] = new float[n];
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (j == 0) {
               T[i][j] = r[i][j]; 
            } else {
                T[i][j] = r[i][j] + T[i][j-1];
            }
        }
    }
    
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {
    //         cout << T[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------" << endl;

    int nR = B/2+1;
    
    float **R = new float*[m];
    for (int i = 0; i < m; i++) {
        R[i] = new float[nR];
    }
    int **S = new int*[m];
    for (int i = 0; i < m; i++) {
        S[i] = new int[nR];
    }
    
    for (int b = 0; b < nR; b++) {
        //cout << b << endl;
        if (b < n) {
            R[0][b] = T[0][b];
            S[0][b] = b;
        } else {
            R[0][b] = T[0][n-1];
            S[0][b] = n-1;
        }
    }
    
    for (int i = 1; i < m; i++) {
        for (int b = 0; b < nR; b++) {
            int max_v = -100;
            int arg_v = -1;
            for (int j = 0; j < n; j++) {
                int idx = b-j;
                if (idx >= 0) {
                    int v = R[i-1][idx] + T[i][j];
                    if (v > max_v) {
                        max_v = v;
                        arg_v = j;
                    }
                }
            }
            
            R[i][b] = max_v;
            S[i][b] = arg_v;
        }
    }
    
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < nR; j++) {
    //         cout << R[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------" << endl;
    
    int max = R[m-1][nR-1];
    int cost = 0;
    int j = nR-1;
    int *js = new int[m];
    for (int i = m-1; i >= 0; i--) {
        int m_j = S[i][j];
        cost += 2*m_j;
        js[i] = m_j;
        j = j - S[i][j];
    }
    
    cout << "max=" << max << ", cost=" << cost;
    cout << ", js=";
    for (int i = 0; i < m; i++) {
        cout << js[i] << ",";
    }
    cout << endl;
}

void run_oprsc_v1(string input, int B, int m, int n) {
    float **r = new float*[m];
    for (int i = 0; i < m; i++) {
        r[i] = new float[n];
    }
    
    string filename = "input/input-" + input + ".csv";
    ifstream data(filename);
    string line;
    int i = 0;
    while(getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        int j = 0;
        while(getline(lineStream, cell, ',')) {
            r[i][j++] = stoi(cell);
        }
        i++;
    }
    
    
    float **T = new float*[m];
    for (int i = 0; i < m; i++) {
        T[i] = new float[n];
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (j == 0) {
               T[i][j] = r[i][j]; 
            } else {
                T[i][j] = r[i][j] + T[i][j-1];
            }
        }
    }
    
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {
    //         cout << T[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------" << endl;

    int nR = B/2+1;
    
    float **R = new float*[m];
    for (int i = 0; i < m; i++) {
        R[i] = new float[nR];
    }
    int **S = new int*[m];
    for (int i = 0; i < m; i++) {
        S[i] = new int[nR];
    }
    
    for (int b = 0; b < nR; b++) {
        //cout << b << endl;
        if (b < n) {
            R[0][b] = T[0][b];
            S[0][b] = b;
        } else {
            R[0][b] = T[0][n-1];
            S[0][b] = n-1;
        }
    }
    
    for (int i = 1; i < m; i++) {
        // cout << i << endl;
        
        for (int b = 0; b < nR; b++) {
            if (b < i) {
                R[i][b] = -1;
                S[i][b] = -1;
            } else if (b == i) {
                R[i][b] = 0;
                S[i][b] = 0;
            } else {
                int max_v = -100;
                int arg_v = -1;
                for (int j = 0; j < n; j++) {
                    int idx = b-j-1;
                    if (idx >= 0) {
                        if (R[i-1][idx] != -1) {
                            int v = R[i-1][idx] + T[i][j];
                            if (v > max_v) {
                                max_v = v;
                                arg_v = j;
                            }
                        }   
                    }
                }
                
                R[i][b] = max_v;
                S[i][b] = arg_v;
            }
        }
    }
    
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < nR; j++) {
    //         cout << R[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------" << endl;
    
    int max = -1;
    int max_m = -1;
    for (int i = 0; i < m; i++) {
        if (R[i][nR-1] > max) {
            max = R[i][nR-1];
            max_m = i;
        }
    }
    
    int cost = 2*max_m;
    int j = nR-1;
    for (int i = max_m; i >= 0; i--) {
        int m_j = S[i][j];
        cost += 2*m_j;
        j = j - (S[i][j]+1);
    }
    
    cout << "max=" << max << ", cost=" << cost << endl;
}
