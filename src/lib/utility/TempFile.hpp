#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

class TempFile {
public:
    TempFile(std::string tmpdir = "", const std::string& contents = "") {
        if (tmpdir.empty()) {
            if ( (tmpdir = getenv("TMPDIR")).empty() )
                tmpdir = "/tmp";
        }
        _path = (boost::filesystem::path(tmpdir) / "BassovacTemp.XXXXXX").string();
        int fd = mkstemp(&_path[0]);
        if (fd < 0)
            throw std::runtime_error("Failed to create temporary file");
        close(fd);
        if (!contents.empty()) {
            std::ofstream f(_path.c_str());
            f << contents;
            f.close();
        }
    }

    ~TempFile() {
        boost::filesystem::remove(_path);
    }

    const std::string& path() const {
        return _path;
    }

protected:
    std::string _path;
};

class TempDir {
public:
    TempDir(const std::string& tmpl = "/tmp/BassovacTempdir.XXXXXX")
        : _path(tmpl)
    {
        using boost::format;
        if (mkdtemp(&_path[0]) != &_path[0])
            throw std::runtime_error(str(format("Failed to create temporary directory '%1%': %2%") %_path %strerror(errno)));
    }

    ~TempDir() {
        boost::filesystem::remove_all(_path);
    }

    std::unique_ptr<TempFile> tempFile(const std::string& contents = "") {
        return std::unique_ptr<TempFile>(new TempFile(path(), contents));
    }

    std::string subpath(std::string const& leaf) const {
        return (boost::filesystem::path(_path) / leaf).string();
    }

    const std::string& path() const {
        return _path;
    }

protected:
    std::string _path;
};
