/**
 * \file ProlificIOHelpers.hpp
 * \author Ted Holzman
 * \date Apr 11, 2012
 * \par Library:
 * Part of the galosh prolific library.
 * \brief contains classes to help in the scanning and parsing of profiles
 */

#ifndef PROLIFIC_IO_HELPERS_HPP_
#define PROLIFIC_IO_HELPERS_HPP_

#include <boost/iostreams/char_traits.hpp> // EOF, WOULD_BLOCK
#include <boost/iostreams/concepts.hpp>    // input_filter
#include <boost/iostreams/operations.hpp>  // get
#include <string>

namespace io = boost::iostreams;
using namespace io;

namespace galosh {

class input_comment_diversion_filter : public input_filter
{
public:
    explicit input_comment_diversion_filter(char comment_char = '#')
      : comment_char_(comment_char), skip_(false)
       { raw_comments_ = std::string("");}
    /// \todo very primitive.  What if it runs into a '#' in quotes?  Must make this into (i.e. find/use) a real scanner some day.
    template<typename Source>
    int get (Source& src)
    {
        int c;
        while (true) {
            if ((c = io::get(src)) == EOF || c == WOULD_BLOCK) break;
            skip_ =
               c == comment_char_ ?
                 true :
                    c == '\n' ?
                       false :
                       skip_;
            if (!skip_)
                break;
            if(c == comment_char_) {
              raw_comments_.append(" ");
            } else {
              raw_comments_.append(const_cast<char *>((char *)(&c)),1L);
            }
        }
        return c;
    }

    std::string get_raw_comments ()
    {
      return raw_comments_;
    } // get_raw_comments

    template<typename Source>
    void close(Source&) { skip_ = false; raw_comments_ = std::string("");} // close;
private:
    char comment_char_;
    bool skip_;
    std::string raw_comments_;
}; // of class input_comment_diversion_filter

}

#endif /* PROLIFIC_IO_HELPERS_HPP_ */
