#pragma once
// Reusable 2-bit packed-sequence primitives (shared by fqdup/dart/fastst).
// Encoding: A=0(00) C=1(01) G=2(10) T=3(11), 4 bases/byte, MSB-first —
// identical layout to fqdup src/derep_detail/packed_ops.hpp (extract_packed_part,
// packed_hamming). GC bases (C,G) are exactly those where the two bits differ.

#include <cstdint>

namespace taph {

// GC count of one packed byte (4 bases). Per lane the low bit of (x ^ (x>>1))
// is b0^b1, set iff the base is C(01) or G(10). Masking 0x55 keeps the 4 low
// lane bits; popcount gives the GC total for the byte.
inline uint8_t gc_count_byte(uint8_t x) {
    return static_cast<uint8_t>(__builtin_popcount((x ^ (x >> 1)) & 0x55u));
}

// Count of C+G over packed 2-bit bases in [start, start+len).
// Fast path: whole aligned bytes via gc_count_byte. Unaligned head/tail fall
// back to per-base, exactly like packed_hamming's shift handling.
inline uint16_t gc_count_packed(const uint8_t* packed, int start, int len) {
    uint16_t gc = 0;
    for (int i = 0; i < len; ) {
        int pos    = start + i;
        int byte_i = pos >> 2;
        int bit_i  = pos & 3;
        if (bit_i == 0 && i + 4 <= len) {
            gc = static_cast<uint16_t>(gc + gc_count_byte(packed[byte_i]));
            i += 4;
        } else {
            uint8_t sh   = static_cast<uint8_t>(6 - 2 * bit_i);
            uint8_t base = (packed[byte_i] >> sh) & 0x3u;
            gc = static_cast<uint16_t>(gc + ((base ^ (base >> 1)) & 1u));
            ++i;
        }
    }
    return gc;
}

}  // namespace taph
