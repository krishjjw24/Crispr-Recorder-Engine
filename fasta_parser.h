#ifndef FASTA_PARSER_HPP
#define FASTA_PARSER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <cctype>

namespace crispr {

/**
 * Loads a single- or multi-record FASTA file into one contiguous uppercase string.
 * Strips '>' header lines and whitespace/newlines for fast sliding-window scans.
 */
inline std::string load_fasta(const std::string& filepath) {
    std::ifstream infile(filepath, std::ios::in | std::ios::binary);
    if (!infile.is_open()) {
        throw std::runtime_error("Error: Could not open FASTA file: " + filepath);
    }

    // Pre-allocate memory based on file size to avoid repeated vector/string reallocations
    infile.seekg(0, std::ios::end);
    std::streamsize file_size = infile.tellg();
    infile.seekg(0, std::ios::beg);

    std::string sequence;
    if (file_size > 0) {
        sequence.reserve(static_cast<size_t>(file_size));
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '>') {
            continue; // Skip FASTA header lines
        }
        for (char c : line) {
            char upper_c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (upper_c == 'A' || upper_c == 'C' || upper_c == 'G' || upper_c == 'T' || upper_c == 'N') {
                sequence.push_back(upper_c);
            }
        }
    }
    return sequence;
}

/**
 * Returns the Watson-Crick complement of a single nucleotide.
 */
constexpr char complement_base(char c) noexcept {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'N';
    }
}

/**
 * Generates the reverse complement of a DNA sequence (5' -> 3')
 * so we can scan the reverse (-) strand in a single forward pass.
 */
inline std::string reverse_complement(std::string_view dna) {
    std::string rev_comp;
    rev_comp.resize(dna.size());

    size_t n = dna.size();
    for (size_t i = 0; i < n; ++i) {
        rev_comp[i] = complement_base(dna[n - 1 - i]);
    }
    return rev_comp;
}

} // namespace crispr

#endif // FASTA_PARSER_HPP