#include <iostream>
#include <fstream>
#include <cstdio>

#include <seqan/stream.h>

// This template function reads the contents from the given Stream in and
// writes it out to std::cout

template <typename TStream>
void doReading(TStream & stream)
{
    typename seqan::DirectionIterator<TStream, seqan::Input>::Type reader;
    seqan::CharString buffer;
    reserve(buffer, 1000);

    reader = directionIterator(stream, seqan::Input());

    while (!atEnd(reader))
    {
        clear(buffer);
        read(buffer, reader, capacity(buffer));
        write(std::cout, buffer);
    }
}

// This template function writes out "Hello World!\n" to the given Stream.

template <typename TStream>
void doWriting(TStream & stream)
{
    seqan::CharString buffer = "Hello World!\n";
    write(stream, buffer);
}

// The main function parses the command line, opens the files in the
// appropriate modes with the appropriate stream types and then calls either
// doWriting() or doReading().

int main(int argc, char const ** argv)
{
    if (argc != 3)
    {
        std::cerr << "USAGE: " << argv[0] << " [r|w] FILENAME\n";
        return 1;
    }

    // Check second argument.
    if (seqan::CharString(argv[1]) != "r" && seqan::CharString(argv[1]) != "w")
    {
        std::cerr << "ERROR: " << argv[1] << " is not a valid operation name.\n";
        return 1;
    }
    bool doRead = (seqan::CharString(argv[1]) == "r");

    // Branches for stream and operation type.
    if (doRead)  // reading
    {
        seqan::VirtualStream<char, seqan::Input> stream;
        if (!open(stream, argv[2]))
        {
            std::cerr << "ERROR: Could not open " << argv[2] << "\n";
            return 1;
        }
        doReading(stream);
    }
    else
    {
        seqan::VirtualStream<char, seqan::Output> stream;
        if (!open(stream, argv[2]))
        {
            std::cerr << "ERROR: Could not open " << argv[2] << "\n";
            return 1;
        }
        doWriting(stream);
    }

    return 0;
}
