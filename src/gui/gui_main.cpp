#include <iostream>
#include <QtGui>
#include <gui_mainWindow.h>

int main(int argc, char * argv[])
{
  typedef std::ostream Log;
  Log & log = std::cout;

  QApplication app(argc, argv);

  MainWindow mainWindow;
  mainWindow.show();
  mainWindow.raise();
  mainWindow.activateWindow();

  //QString fileName = QFileDialog::getOpenFileName(0,"Open", ".");

  //std::cout << "Filename: " << fileName.toStdString() << std::endl;

  return app.exec();
}
