#include "bvprob/FastaReader.hpp"
#include "utility/TempFile.hpp"

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

using namespace std;

class TestFastaReader : public testing::Test {
public:
    TempDir tmpdir;
};

TEST_F(TestFastaReader, noData) {

    // file not found
    ASSERT_THROW(FastaReader("/_____/no_fasta/files_here.:p"), runtime_error);

    string data = ">1\nAA\n";
    auto fasta = tmpdir.tempFile(data);
    FastaReader reader(fasta->path());

    // sequence not found
    ASSERT_THROW(reader.sequence("3", 0), length_error);
}

TEST_F(TestFastaReader, multiSequence) {
    string data =
        ">1\n"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
        ">2\n"
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
        ">3\n"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"
        ">4\n"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";

    auto fasta = tmpdir.tempFile(data);
    FastaReader reader(fasta->path());

    // read sequentially
    for (unsigned i = 0; i < 40; ++i)
        ASSERT_EQ('A', reader.sequence("1", i));
    for (unsigned i = 0; i < 40; ++i)
        ASSERT_EQ('C', reader.sequence("2", i));
    for (unsigned i = 0; i < 40; ++i)
        ASSERT_EQ('G', reader.sequence("3", i));
    for (unsigned i = 0; i < 40; ++i)
        ASSERT_EQ('T', reader.sequence("4", i));

    // interleave chromosome access
    for (unsigned i = 0; i < 40; ++i) {
        ASSERT_EQ('A', reader.sequence("1", i));
        ASSERT_EQ('C', reader.sequence("2", i));
        ASSERT_EQ('G', reader.sequence("3", i));
        ASSERT_EQ('T', reader.sequence("4", i));
    }
}

TEST_F(TestFastaReader, offset) {
    string data = ">1\n" "ACGT\n";
    auto fasta = tmpdir.tempFile(data);
    FastaReader reader(fasta->path());
    ASSERT_EQ('A', reader.sequence("1", 0));
    ASSERT_EQ('C', reader.sequence("1", 1));
    ASSERT_EQ('G', reader.sequence("1", 2));
    ASSERT_EQ('T', reader.sequence("1", 3));
}

TEST_F(TestFastaReader, readPastEndOfChromosome) {
    string data = ">1\n" "ACGT\n";
    auto fasta = tmpdir.tempFile(data);
    FastaReader reader(fasta->path());
    ASSERT_THROW(reader.sequence("1", 4), length_error);
}
