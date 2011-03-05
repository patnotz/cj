#ifndef GUI_MAINWINDOW_H_
#define GUI_MAINWINDOW_H_

#include <QMainWindow>
class QAction;
class QLabel;
class QLineEdit;
class QPushButton;
class QTextBrowser;

#include <json/json.h>

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();

protected:
	void closeEvent(QCloseEvent * event);

private slots:
  	void openInputMesh();
	void openOutputMesh();
	void run();

private:
	void createWidgets();
	void createActions();
	void createMenus();
	void createContextMenu();
	void createToolBars();
	void createStatusBar();

	QWidget *centralWidget;

	QLabel *inputMeshLabel;
	QLineEdit *inputMeshLineEdit;
	QPushButton *inputMeshOpenPushButton;

	QLabel *outputMeshLabel;
	QLineEdit *outputMeshLineEdit;
	QPushButton *outputMeshOpenPushButton;

	QTextBrowser *outputTextBrowser;

	QPushButton *runPushButton;
	QAction * quitAction;

	Json::Value config;
};

#endif /* GUI_MAINWINDOW_H_ */
