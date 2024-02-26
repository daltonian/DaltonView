#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include <cmath>
#include <thread>
using namespace std;


// Create function for reading lines from output file
vector<string> readLinesFromFile(const string& filePath) {
    vector<string> lines;
    ifstream file(filePath);

    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
    } else {
        cerr << "Unable to open file: " << filePath << endl;
    }

    return lines;
}


// Create function for finding an index for a given Regex pattern
vector<int> findIndexUsingRegex(const vector<string>& lines, const regex& pattern, int start = 0, int end = -1) {
    if (end == -1) {
        end = lines.size();
    }

    vector<int> indices;
    for (size_t i = start; i < end; ++i) {
        if (regex_search(lines[i], pattern)) {
            indices.push_back(static_cast<int>(i));
        }
    }
    return indices;
}

// Use index finder to get indices for data values (there are many per file)
void getDataIndices(const vector<string>& lines, vector<int>* eV_index, vector<int>* OscStr_index) {
    regex p ("TDDFT Excitation Energies"); // Find index of where output data begins
    int start = findIndexUsingRegex(lines,p)[0];

    p = ("SETman timing");  // Set the pattern, p, for regex test
    int end = findIndexUsingRegex(lines,p)[0]; // Find index of where output data ends

    p = ("Excited state"); // Set the pattern, p, for regex test
    *eV_index = findIndexUsingRegex(lines,p,start,end);  // Find indices of excited-state energy values

    p = ("Strength");  // Set the pattern, p, for regex test
    *OscStr_index = findIndexUsingRegex(lines,p,start,end);  // Find indices of oscillator strength values
}

// Fetch numerical values from indices
void getDataValues(const vector<string>& lines, vector<int> eV_index, vector<int> OscStr_index,
 vector<double>* eV_values, vector<double>* wl_nm_values, vector<double>* OscStr_values) {
    regex pattern ("\\b\\d*\\.\\d+\\b");
    smatch match;
    for (size_t i=0; i < eV_index.size(); i++){
        regex_search(lines[eV_index[i]],match,pattern);
        eV_values->push_back(stod(match[0]));
        wl_nm_values->push_back((1239.8) / stod(match[0]));
    }

    for (size_t i=0; i < OscStr_index.size(); i++){
        regex_search(lines[OscStr_index[i]],match,pattern);
        OscStr_values->push_back(stod(match[0]));
    }
}

// Create new struct so it can be returned from function
struct Spectrum {
        vector<double> wavelength;
        vector<double> absorbance;
};

// Apply fit functions to sticks
Spectrum createAbsorptionSpectrum(vector<double> wl_nm_values, vector<double> OscStr_values) {
    vector<double> wl_range;
    for (float i=100; i<700;i=i+0.5) {
        wl_range.push_back(i);
    }

    cout << "Enter the hwhm of your peaks in nm: " << endl;
    float hwhm;
    cin >> hwhm;
    
    double pi = 3.14159265358979323846;

    Spectrum result;
    result.wavelength = wl_range;

    for (double x : wl_range) {
        double sum = 0.0;
        for (size_t i=0;i<wl_nm_values.size();i++) {
            if (wl_nm_values[i] > 200) {
                double nm = wl_nm_values[i];
                double OS = OscStr_values[i];
                sum += OS * exp(-0.5 * pow((x - nm), 2) / (pow(hwhm, 2))) / (hwhm * sqrt(2 * M_PI));
            }
        }
        result.absorbance.push_back(sum);
    }
    return result;
}

void exitProgram() {
    cout << "**********************************************" << endl;
    cout << "*       Thank you for using DaltonView       *" << endl;
    cout << "**********************************************" << endl;
    
    cout << "Program will close in 5 seconds." << endl;
    chrono::seconds duration(5);
    this_thread::sleep_for(duration);
}


void printHeader() {
    cout << "**********************************************" << endl;
    cout << "*          DaltonView, CMD Edition           *" << endl;
    cout << "**********************************************" << endl;
    cout << "* This version is optimized for Q-Chem files *" << endl;
    cout << " " << endl;
}

string getFilePath() {
    string filePath;
    cout << "What is the path for the TDDFT file?" << endl;
    getline(std::cin, filePath);
    cout << "Read input path as " << filePath << endl;
    // Remove extra quotation marks, if necessary
    filePath.erase(remove(filePath.begin(), filePath.end(), '\"'), filePath.end());
    return filePath;
}

int main() {

    printHeader();

    // Take input from user for file path
    string filePath = getFilePath();

    
    // Intake lines into file as string vector
    vector<string> lines = readLinesFromFile(filePath);
    
    // Exit if nothing actually loaded
    if (lines.empty()) {
        std::cout << "No file could be loaded." << std::endl;
        exitProgram();
        return 1;
    }

    vector<int> eV_index;
    vector<int> OscStr_index;

    // Find indicies in 'lines' of data values
    // Pointers are passed to avoid needing to return vectors
    getDataIndices(lines,&eV_index,&OscStr_index);

    vector<double> eV_values;
    vector<double> wl_nm_values;
    vector<double> OscStr_values;

    // Use data indicies to fetch data values with regex
    // Pointers used for data value vectors 
    getDataValues(lines,eV_index,OscStr_index,&eV_values,&wl_nm_values,&OscStr_values);

    Spectrum result = createAbsorptionSpectrum(wl_nm_values,OscStr_values);

    string outputPath = filePath.replace(filePath.length() - 4, 4, "_spectrum.txt");
    ofstream outputFile(outputPath);

    if (outputFile.is_open()) {
        // Write headers
        outputFile << "wavelength\tabsorbance\n";

        // Write data to the file with tabs
        for (size_t i = 0; i < result.wavelength.size(); ++i) {
            outputFile << result.wavelength[i] << "\t" << result.absorbance[i] << "\n";
        }

        // Close the file
        outputFile.close();
        cout << "Data exported to " << outputPath << endl;
    } else {
        cerr << "Unable to open file for writing\n";
    }
    exitProgram();
    return 0;
}
