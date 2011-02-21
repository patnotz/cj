#include <QtGui>
#include <gui_mainWindow.h>

MainWindow::MainWindow() {
	outputTextBrowser = new QTextBrowser;
	setCentralWidget(outputTextBrowser);

	outputTextBrowser->append("Hello, World!");

	createActions();
	createMenus();
	createContextMenu();
	createToolBars();
	createStatusBar();
}

void MainWindow::createActions() {
	quitAction = new QAction("&Quit", this);
	quitAction->setShortcut(QKeySequence("Ctrl+Q"));
	quitAction->setStatusTip("Quit the application that is not called Cj");
	connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));
}

void MainWindow::createMenus() {

}

void MainWindow::createContextMenu() {

}

void MainWindow::createToolBars() {

}

void MainWindow::createStatusBar() {

}

void MainWindow::closeEvent(QCloseEvent * event) {

}

void MainWindow::openInputMesh() {

}

void MainWindow::openOutputMesh() {

}

void MainWindow::run() {

}
