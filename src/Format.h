#ifndef FORMAT_H
#define FORMAT_H

static inline void GetFormatSubstrings(std::vector<std::string>& subStrings, std::string templateStr)
{
    std::string delimiter = "{}";
    std::string templateStrCopy = templateStr;
    size_t pos = 0;
    std::string token;
    while ((pos = templateStrCopy.find(delimiter)) != std::string::npos)
    {
        token = templateStrCopy.substr(0, pos);
        subStrings.push_back(token);
        templateStrCopy.erase(0, pos + delimiter.length());
    }
    subStrings.push_back(templateStrCopy);
}

static inline std::vector<std::string> StringSplit(std::string templateStr, std::string delimiter)
{
    std::vector<std::string> subStrings;
    std::string templateStrCopy = templateStr;
    size_t pos = 0;
    std::string token;
    while ((pos = templateStrCopy.find(delimiter)) != std::string::npos)
    {
        token = templateStrCopy.substr(0, pos);
        subStrings.push_back(token);
        templateStrCopy.erase(0, pos + delimiter.length());
    }
    subStrings.push_back(templateStrCopy);
    return subStrings;
}

template <typename T> static inline void strformat_recursive (std::string& templateStr, std::ostringstream& strstream, std::vector<std::string>& subStrings, int& lev, T t)
{
    strstream << subStrings[lev] << t;
    lev++;
}

template <typename T, typename... Ts> static inline void strformat_recursive (std::string& templateStr, std::ostringstream& strstream, std::vector<std::string>& subStrings, int& lev, T t, Ts... ts)
{
    strstream << subStrings[lev] << t;
    lev++;
    strformat_recursive(templateStr, strstream, subStrings, lev, ts...);
}

template <typename... Ts> static inline std::string strformat (std::string templateStr, Ts... ts)
{
    std::ostringstream strstream;
    std::vector<std::string> subStrings;
    GetFormatSubstrings(subStrings, templateStr);
    if ((sizeof...(Ts))!=(subStrings.size()-1))
    {
        CallError("Formatted string \"" + templateStr + "\" has wrong number of placeholders!");
    }
    int lev = 0;
    strformat_recursive(templateStr, strstream, subStrings, lev, ts...);
    strstream << subStrings[lev];
    return strstream.str();
}
#endif