#include <iostream>
#include <string>
#include <string_view>
#include <cassert>
#include "../include/fasta_parser.h"

int main() {
    try {
        // 1. Unit-test Reverse Complement (Rosalind REVC verification)
        std::string_view sample = "AAAACCCGGT";
        std::string revc = crispr::reverse_complement(sample);
        assert(revc == "ACCGGGTTTT");
        std::cout << "[PASS] Unit Test 1: Reverse complement of " << sample
                  << " -> " << revc << "\n";

        // 2. Load the FASTA genome file
        const std::string filepath = "C:/Users/User/Downloads/CPP/crispr-recorder-engine/data/test_ecoli.fasta";
        std::string forward_strand = crispr::load_fasta(filepath);
        std::string reverse_strand = crispr::reverse_complement(forward_strand);

        std::cout << "[PASS] Unit Test 2: Loaded FASTA sequence ("
                  << forward_strand.size() << " bp)\n";
        std::cout << "       5' Forward Preview: "
                  << std::string_view(forward_strand).substr(0, 40) << "...\n";
        std::cout << "       5' Reverse Preview: "
                  << std::string_view(reverse_strand).substr(0, 40) << "...\n\n";

        // 3. Zero-copy std::string_view 23-nt sliding window preview (20-nt spacer + 3-nt NGG PAM)
        std::string_view fwd_view(forward_strand);
        size_t ngg_count_fwd = 0;

        std::cout << "--- Scanning Forward Strand for 3' NGG (SpCas9) Windows ---\n";
        for (size_t i = 0; i + 23 <= fwd_view.size(); ++i) {
            // Check if positions [i+21, i+22] are 'GG' (meaning [i+20, i+22] is 'NGG')
            if (fwd_view[i + 21] == 'G' && fwd_view[i + 22] == 'G') {
                std::string_view spacer = fwd_view.substr(i, 20);   // Zero-copy 20-nt slice!
                std::string_view pam    = fwd_view.substr(i + 20, 3);
                ++ngg_count_fwd;

                if (ngg_count_fwd <= 5) {
                    std::cout << "  Hit #" << ngg_count_fwd
                              << " | Pos: " << i
                              << " | Spacer (20-nt): " << spacer
                              << " | PAM: " << pam << "\n";
                }
            }
        }
        std::cout << "Total Forward-Strand NGG Candidates Found: " << ngg_count_fwd << "\n";

    } catch (const std::exception& ex) {
        std::cerr << ex.what() << "\n";
        return 1;
    }
    return 0;
}