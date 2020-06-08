class ProgressBar
{
public:

  ProgressBar(unsigned int max, std::string message = std::string("")) :
    _max(max), _progress(0), _counter(0), _message(message),
    _bar(50, ' ')
  {
    _bar[0] = '#';
  }

  ~ProgressBar()
  {
    std::cout << std::endl;
  }

  void progress()
  {
    unsigned int currentProgress = (_counter*100)/(_max - 1);

    if (_progress < currentProgress)
    {
      _progress = currentProgress;

      std::fill_n(_bar.begin(), _progress/2, '#');

      std::cout << "\r" << _message << "  [" << _bar << "] "
                << _progress << "%" << std::flush;
    }

    ++_counter;
  }

  void operator++()
  {
    progress();
  }

  void operator++(int)
  {
    progress();
  }

protected:

  unsigned int _max;
  unsigned int _progress;
  unsigned int _counter;

  std::string _message;
  std::string _bar;
};
