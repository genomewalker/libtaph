#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <string>

namespace taph {

// Sequential reader of read sequences. The host tool (fqdup, dart) supplies it, so
// libtaph needs no FASTQ or compression code.
class SeqReader {
public:
    virtual ~SeqReader() = default;
    // Next sequence in file order; false at end of input.
    virtual bool next(std::string& seq) = 0;
};

// Opens `path` at its first record. `threads` is a decompression hint.
using SeqReaderOpener = std::function<std::unique_ptr<SeqReader>(const std::string& path, size_t threads)>;

} // namespace taph
