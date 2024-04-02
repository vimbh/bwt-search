// Written by Vimukthi Herath for COMP9319 23T2

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_set>

bool MSB_isOne(char &value) {
    // 1 << 7 -> 10000000
    return value & (1 << 7);
}

std::pair<std::vector<int64_t>, std::vector<std::unordered_map<uint8_t, int>>> createIndex(const std::string& T, int& n) {

    uint8_t max_char = 127;

    // Create C Table (This basically represents F in FL mapping).

    std::vector<int64_t> C_table(max_char, 0);

    for (auto& ch : T) {
        C_table[static_cast<uint8_t>(ch)]++;
    }

    int64_t sum = 0;
    for (uint8_t i = 0; i < max_char; ++i) {
        sum += C_table[i];
        C_table[i] = sum - C_table[i];
    }

    // Create Occurrence Table: Occ[index][char].

    int intervals = static_cast<int>(T.length()) / n;
    std::vector<std::unordered_map<uint8_t, int>> occ_table(intervals);

    // Init all characters in each map
    for (size_t i = 0; i < T.length(); ++i) {
        if (occ_table[0].find(static_cast<uint8_t>(T[i])) == occ_table[0].end()) {
            occ_table[0][static_cast<uint8_t>(T[i])] = 0;
        }
    }
    
    std::unordered_map<uint8_t, int> temp_map;
    if (!occ_table.empty() && !occ_table[0].empty()) {
        temp_map = occ_table[0];
    }
    temp_map[static_cast<uint8_t>(T[0])] = 1;
    // Loop over for remaining prefixes
    for (size_t i = 1; i < T.length(); ++i) {
        temp_map[static_cast<uint8_t>(T[i])] += 1;
        if ((i + 1) % n == 0) {
            occ_table[(i + 1) / n - 1] = temp_map;
        }
    }


    return std::make_pair(C_table, occ_table);
}

int getOccurrences(char c, const std::vector<std::unordered_map<uint8_t, int>>& occ, int row, const std::string& T, int n) {

    int adj_Row = (row + 1) / n - 1;

    if (adj_Row < 0) {
        adj_Row = 0;
    }

    int closest_Row = (adj_Row + 1) * n - 1;
    int dist_Row = std::abs(row - closest_Row);
    int req_occurrences = occ[adj_Row].at(static_cast<uint8_t>(c));

    char ch;

    // Closest row is before the needed one, increment count as required.
    if (closest_Row < row) {

        for (int i = 1; i <= dist_Row; i++) {
            ch = T[closest_Row + i];
            if (ch == c) {
                req_occurrences += 1;
            }
        }
    
    }

    // Closest row is after the needed one, decrement count as required.
    if (closest_Row > row) {

        // Check initial row
        if (T[closest_Row] == c) {
            req_occurrences -= 1;
        }
            // pretty sure correct to do <, not <= (check again later)
        for (int i = 1; i < dist_Row; i++) {
            ch = T[closest_Row - i];
            if (ch == c) {
                req_occurrences -= 1;
            }
        }
    
    }

    return req_occurrences;
}

int getRow(int i, char c, const std::vector<int64_t>& C, const std::vector<std::unordered_map<uint8_t, int>>& occ, const std::string& T, int n) {

    int req_i = i;
    // If i < 0, occ[c-1] doesn't exist; we don't need to adjust for previous occurrences of c.
    if (req_i >= 0) {
        
        if ((req_i + 1) % n != 0) {
            int occ_req = getOccurrences(c, occ, req_i, T, n);
            req_i = C[static_cast<uint8_t>(c)] + occ_req;
        } else {
            req_i = C[static_cast<uint8_t>(c)] + occ[(req_i + 1) / n - 1].at(static_cast<uint8_t>(c));
        }
    } else {
        req_i = C[static_cast<uint8_t>(c)];
    }

    return req_i;
}



int main(int argc, char* argv[]) {

    // Use (var & MSB) to flip MSB of value to 0; defined once to avoid redundancy
    uint8_t msb = 0; // 0000
    msb -= 1; // 1111
    msb >>= 1; // 0111


    if (argc < 3) {
        std::cerr << "Error: Please provide 3 arguments." << std::endl;
        return 1; 
    }

    std::string path_to_text = argv[1];
    std::string path_to_index = argv[2];
    std::ifstream inFile(path_to_text, std::ios::binary);

    // Required Pattern
    std::string P = argv[3];

    //BWT Transformed string
    std::string T = "";

    // Load string T into memory
    if(inFile.is_open()) {
        char c; 
        char prev_char;
        while(inFile.get(c)) {
            if (!MSB_isOne(c)) {
                T += c;
                prev_char = c;
            } else {
                uint8_t prev = c & msb;
                auto curr = prev;
                char nextVal = inFile.peek();
                
                while (MSB_isOne(nextVal) && nextVal != EOF) {
                    inFile.get(c);
                    uint8_t next = c & msb;
                    curr = (next << 8) | prev;
                    prev = next;
                    nextVal = inFile.peek();
                }
                auto currVal = static_cast<unsigned>(curr);
                for (unsigned i = 0; i < currVal+2; i++) {
                    T += prev_char;
                }                
            }

        }

    } else {
        std::cerr << "Error: Unable to open the file." << std::endl;
        return 1;
    }
    inFile.close();


    // Select some interval n for occurance table
    // Scale at 10% of text
    int n = 0.1 * T.length();

    // Create Index
    const auto [C, occ] = createIndex(T, n);

    // Handle the case when occ vector is empty
    if (occ.empty()) {
        return 0;
    }

    // Backwards search from last char in P.
    auto i = P.size() - 1;
    char c = P[i];
    std::vector<std::string> results;

    // If the char is not in the text, no match
    if (occ[0].find(c) == occ[0].end()) {
        return 0;
    }

    // Initial range for P[0]
    int First = C[static_cast<uint8_t>(c)];
    int Last = C[static_cast<uint8_t>(c) + 1] - 1;

    // Loop over chars in P
    while (i > 0) {

        c = P[i - 1];

        if (occ[0].find(c) == occ[0].end()) {
            return 0;
        }

        First = getRow(First-1,c,C,occ,T,n);
        Last = getRow(Last,c,C,occ,T,n) - 1;

        // If Last-First returns lt. 0, there is no correct FL mapping; no match.
        if (Last - First < 0) {
            return 0;
        }

        i = i - 1;

    }

    // Init P2, equal to nextID.
    std::string P2 = "";
    int i_at_match;
    int j;
    
    // Each character in the range [First,Last] corresponds to a match for P.
    for (int i = First; i <= Last; ++i) {

        i_at_match = i;
        j = i;

        // matchFront represents matched string from the opening bracket of the relevant ID '[' to the last matched char in P; P[0].
        // Run backward search from P[0] to '[' to retrieve matchFront.
        std::string matchFront(1, T[j]);

        char c = T[j];

        while (c != '[') {
            j = getRow(j-1, c, C, occ, T, n);

            matchFront = T[j] + matchFront;
            c = T[j];

        } 

        // Retrieve current ID and following ID
        std::string id;
        for (char ch : matchFront.substr(1)) {
            if (ch != ']') {
                id += ch;
            } else {
                break;
            }
        }
    

        std::string next_id = "[" + std::to_string(std::stoi(id) + 1) + "]";
        id = "[" + id + "]";

        P2 = next_id;

        auto k = P2.length() - 1;
        c = P2[k];

        // z will be where the next search is from.
        int z;

        // Initial range for P2[0]
        int First2 = C[static_cast<uint8_t>(c)];
        int Last2 = C[static_cast<uint8_t>(c) + 1] - 1;

        bool lastIdInFile = false;
        // ID's may have common suffixes, e.g. [8] & [18]
        // When rows converge (First2 == Last2), we have found the exact match on the next ID.
        while (First2 != Last2) {

            c = P2[k - 1];

            // If we can't retrieve the occurrence of an integer in the ID despite having a match -> next ID doesn't exist.
            if (occ[0].find(static_cast<uint8_t>(c)) == occ[0].end()) {
                lastIdInFile = true;
                break;
            }
            
            First2 = getRow(First2-1,c,C,occ,T,n);
            Last2 = getRow(Last2,c,C,occ,T,n) - 1;

            k = k - 1;

            // If Last2-First2 returns lt. 0, there is no correct FL mapping despite having a match -> next ID doesn't exist.
            if (Last2 - First2 < 0) {
                lastIdInFile = true;
                break;
            }

            // If it's not the last ID in the file, Backward Search starts from First2 = Last.
            z = First2; 

        }

        // Edge case: Matched ID is last in the file.
        // We perform backward search from the first ID in the file to the matched (Last) ID. As the BWT maps rotations,
        // rotating backward from the start of the file means the first ID encountered is the last ID.
        // Need to first determine the first ID in the file, as it could be any non-negative integer.
        if (lastIdInFile) {

            int min = std::numeric_limits<int>::max();
            char c = ']';

            // Initial range for search: All matches on "]"
            int First3 = C[static_cast<uint8_t>(c)];
            int Last3 = C[static_cast<uint8_t>(c) + 1] - 1;

            // Backward Search on all ID's
            int m;
            int t;
            int min_position;
            for (t = First3; t <= Last3; ++t) {

                m = t;
                std::string idString(1,T[m]);
                char c = T[m];
                while (c != '[') {

                    m = getRow(m-1, c, C, occ, T,n);
                    idString = T[m] + idString;
                    c = T[m];
                }

                // Retrieve ID
                std::string curr_id;
                for (char ch : idString.substr(1)) {
                    if (ch != ']') {
                        curr_id += ch;
                    } else {
                        break;
                    }
                }

                // Update smallest ID
                // min_position = row of FM index corresponding to firstID[0]. Next iteration of search finds '[', i.e. the first character in the file.
                if (std::stoi(curr_id) < min) {
                    min = std::stoi(curr_id);
                    min_position = m;

                }

                // Backward Search starts from the first ID until we reach the position (row of FM index) of the last match on the Last ID.
                z = min_position;
            }
            lastIdInFile = false;
        }

        // Once correct ID has been found, continue backward search to resolve our match.
        char c2 = T[z];
        std::string matchEnd = "";

        if (T[z] != '[') {
            matchEnd += T[z];
        }
        
        while (z != i_at_match) {

            z = getRow(z-1, c2, C, occ, T, n);
            if (T[z] != '[' && z != i_at_match) {
                matchEnd = T[z] + matchEnd;
            }
            c2 = T[z];
        }
        std::string fullMatch = matchFront + matchEnd;

        if (std::find(results.begin(), results.end(), fullMatch) == results.end()) {
            results.push_back(fullMatch);
        }
        
        
    }



    std::sort(results.begin(), results.end(), [](const std::string& s1, const std::string& s2) {
        int id1 = std::stoi(s1.substr(1, s1.find(']') - 1));
        int id2 = std::stoi(s2.substr(1, s2.find(']') - 1));
        return id1 < id2;
    });

    // Output Results.
    for (auto& result : results) {
        std::cout << result << std::endl;
    }

    return 0;
}




    