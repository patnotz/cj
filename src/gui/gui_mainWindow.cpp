#include <QtGui>
#include <gui_mainWindow.h>
#include <drive_simulation.h>
#include <log.h>

MainWindow::MainWindow()
{
  createWidgets();
  createActions();
  createMenus();
  createContextMenu();
  createToolBars();
  createStatusBar();
}

void MainWindow::createWidgets()
{
  centralWidget = new QWidget;
  setCentralWidget(centralWidget);

  inputMeshLabel = new QLabel(tr("Input Mesh:"));
  inputMeshLineEdit = new QLineEdit;
  inputMeshOpenPushButton = new QPushButton(tr("Open..."));
  QHBoxLayout *inputMeshHBoxLayout = new QHBoxLayout;
  inputMeshHBoxLayout->addWidget(inputMeshLabel, 0);
  inputMeshHBoxLayout->addWidget(inputMeshLineEdit, 1);
  inputMeshHBoxLayout->addWidget(inputMeshOpenPushButton, 0);
  connect(inputMeshOpenPushButton, SIGNAL(clicked()), this,
    SLOT(openInputMesh()));

  outputMeshLabel = new QLabel("Output Mesh:");
  outputMeshLineEdit = new QLineEdit;
  outputMeshOpenPushButton = new QPushButton(tr("Open..."));
  QHBoxLayout *outputMeshHBoxLayout = new QHBoxLayout;
  outputMeshHBoxLayout->addWidget(outputMeshLabel, 0);
  outputMeshHBoxLayout->addWidget(outputMeshLineEdit, 1);
  outputMeshHBoxLayout->addWidget(outputMeshOpenPushButton, 0);
  connect(outputMeshOpenPushButton, SIGNAL(clicked()), this,
    SLOT(openOutputMesh()));

  outputTextBrowser = new QTextBrowser;
  QHBoxLayout *outputTextBrowserHBoxLayout = new QHBoxLayout;
  outputTextBrowserHBoxLayout->addWidget(outputTextBrowser, 1);

  QHBoxLayout *runHBoxLayout = new QHBoxLayout;
  runPushButton = new QPushButton(tr("Run"));
  runPushButton->setEnabled(true); // FIXME: use some kind of validation
  runHBoxLayout->addWidget(runPushButton, 0, Qt::AlignRight);
  connect(runPushButton, SIGNAL(clicked()), this, SLOT(run()));

  QVBoxLayout *mainWindowVBoxLayout = new QVBoxLayout(centralWidget);
  mainWindowVBoxLayout->addLayout(inputMeshHBoxLayout, 0);
  mainWindowVBoxLayout->addLayout(outputMeshHBoxLayout, 0);
  mainWindowVBoxLayout->addLayout(outputTextBrowserHBoxLayout);
  mainWindowVBoxLayout->addLayout(runHBoxLayout);
}

void MainWindow::createActions()
{
  quitAction = new QAction(tr("&Quit"), this);
  quitAction->setShortcut(QKeySequence("Ctrl+Q"));
  quitAction->setStatusTip(tr("Quit the application that is not called Cj"));
  connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));
}

void MainWindow::createMenus()
{

}

void MainWindow::createContextMenu()
{

}

void MainWindow::createToolBars()
{

}

void MainWindow::createStatusBar()
{

}

void MainWindow::closeEvent(QCloseEvent * event)
{

}

void MainWindow::openInputMesh()
{
  QString fileName = QFileDialog::getOpenFileName(0, tr("Open Input Mesh"),
    QDir::currentPath());
  if (fileName.size() > 0)
  {
    inputMeshLineEdit->setText(fileName);
    // FIXME: config should be updated when the text box
    // changes using signal-slot
    config["mesh-database"] = fileName.toStdString();
  }
}

void MainWindow::openOutputMesh()
{
  QString fileName = QFileDialog::getOpenFileName(0, tr("Open Output Mesh"),
    QDir::currentPath());
  if (fileName.size() > 0)
  {
    outputMeshLineEdit->setText(fileName);
    // FIXME: config should be updated when the text box
    // changes using signal-slot
    config["results-database"] = fileName.toStdString();
  }
}

void MainWindow::run()
{
  // Set the log, error and warning streams. We don't have options for these.
  Log::setLog(std::cout);
  Log::setWarning(std::cout);
  Log::setError(std::cerr);
  drive_simulation(config, 0, NULL);
}
