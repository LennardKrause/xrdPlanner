def main():
    import sys
    from PyQt6 import QtWidgets
    from xrdPlanner.classes import MainWindow
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec())
    
if __name__ == '__main__':
    main()
