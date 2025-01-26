/*
 * XYZ to DXF Converter GUI
 * -------------------------
 * This program provides a graphical user interface (GUI) for the command-line tool `xyz2dxf`.
 * The GUI allows users to specify input parameters and execute the tool without using the terminal.
 *
 * Key Features:
 * - File dialog for selecting input files.
 * - Input fields for required parameters (`minDist`, `precision`, `PDMODE`).
 * - Optional parameters (`gridSpacing`, `maxTPSPoints`) can be specified.
 * - Executes `xyz2dxf` as a child process with the specified parameters.
 *
 * Compilation (Standalone Static Executable):
 * -------------------------------------------
 * Use the following command to compile:
 * g++ -O3 -fopenmp -flto -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o xyz2dxf_gui.exe xyz2dxf_gui.cpp -lkernel32 -lopengl32 -luuid -lcomdlg32
 *
 * Notes:
 * - `-static`: Ensures static linking for a standalone executable.
 * - `-mwindows`: Hides the console window for GUI applications.
 * - `-lcomdlg32`: Links the library required for file dialog functionality.
 *
 * Ensure `xyz2dxf` is in the system PATH for proper execution.
 */

#include <windows.h>
#include <commdlg.h>
#include <string>
#include <sstream>
#include <iomanip>

// Define control IDs for GUI components
#define IDC_INPUT_FILE 1001
#define IDC_MIN_DIST 1002
#define IDC_PRECISION 1003
#define IDC_PDMODE 1004
#define IDC_GRID_SPACING 1005
#define IDC_MAX_TPS_POINTS 1006
#define IDC_BROWSE_BUTTON 1007
#define IDC_RUN_BUTTON 1008

/**
 * openFileDialog:
 * ---------------
 * Opens a file dialog to select an input file.
 *
 * @param hwnd Handle to the parent window.
 * @return A string containing the full path of the selected file, or an empty string if canceled.
 */
std::string openFileDialog(HWND hwnd)
{
    OPENFILENAME ofn; // Structure for file dialog settings
    char szFile[260]; // Buffer to store the selected file path

    ZeroMemory(&ofn, sizeof(ofn));
    ZeroMemory(szFile, sizeof(szFile));

    ofn.lStructSize = sizeof(ofn); // Structure size
    ofn.hwndOwner = hwnd;          // Owner window handle
    ofn.lpstrFile = szFile;        // File buffer
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = "XYZ Files (*.xyz)\0*.xyz\0All Files (*.*)\0*.*\0"; // File type filters
    ofn.nFilterIndex = 1;                                                 // Default filter
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;                    // Validation flags

    // Display the file dialog and return the selected file path
    if (GetOpenFileName(&ofn))
    {
        return std::string(ofn.lpstrFile);
    }
    return ""; // Return an empty string if canceled
}

/**
 * WindowProc:
 * -----------
 * Main window procedure for handling messages sent to the GUI.
 *
 * @param hwnd Handle to the window.
 * @param uMsg The message code.
 * @param wParam Additional message information.
 * @param lParam Additional message information.
 * @return LRESULT Result of message processing.
 */
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    static HWND hInputFile, hMinDist, hPrecision, hPDMODE, hGridSpacing, hMaxTPSPoints;

    switch (uMsg)
    {
    case WM_CREATE:
    {
        // Input File Field
        CreateWindowEx(WS_EX_CLIENTEDGE, "STATIC", "Input File:", WS_VISIBLE | WS_CHILD,
                       20, 20, 200, 20, hwnd, NULL, NULL, NULL);
        hInputFile = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                    150, 20, 250, 25, hwnd, (HMENU)IDC_INPUT_FILE, NULL, NULL);
        CreateWindow("BUTTON", "Browse", WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON,
                     410, 20, 80, 25, hwnd, (HMENU)IDC_BROWSE_BUTTON, NULL, NULL);

        // Minimum Distance Field
        CreateWindowEx(WS_EX_CLIENTEDGE, "STATIC", "Min Dist:", WS_VISIBLE | WS_CHILD,
                       20, 60, 100, 20, hwnd, NULL, NULL, NULL);
        hMinDist = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "5", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                  150, 60, 100, 25, hwnd, (HMENU)IDC_MIN_DIST, NULL, NULL);
        CreateWindow("STATIC",
                     "Minimum distance threshold for filtering out closely spaced points.",
                     WS_VISIBLE | WS_CHILD, 260, 60, 500, 40, hwnd, NULL, NULL, NULL);

        // Precision Field
        CreateWindowEx(WS_EX_CLIENTEDGE, "STATIC", "Precision:", WS_VISIBLE | WS_CHILD,
                       20, 110, 100, 20, hwnd, NULL, NULL, NULL);
        hPrecision = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "2", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                    150, 110, 100, 25, hwnd, (HMENU)IDC_PRECISION, NULL, NULL);
        CreateWindow("STATIC",
                     "Number of decimal places for numerical outputs.",
                     WS_VISIBLE | WS_CHILD, 260, 110, 500, 40, hwnd, NULL, NULL, NULL);

        // PDMODE Field
        CreateWindowEx(WS_EX_CLIENTEDGE, "STATIC", "PDMODE:", WS_VISIBLE | WS_CHILD,
                       20, 160, 100, 20, hwnd, NULL, NULL, NULL);
        hPDMODE = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "3", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                 150, 160, 100, 25, hwnd, (HMENU)IDC_PDMODE, NULL, NULL);
        CreateWindow("STATIC",
                     "Drawing style for points in DXF output (integer code).",
                     WS_VISIBLE | WS_CHILD, 260, 160, 500, 40, hwnd, NULL, NULL, NULL);

        // Grid Spacing Field
        CreateWindowEx(WS_EX_CLIENTEDGE, "STATIC", "Grid Spacing:", WS_VISIBLE | WS_CHILD,
                       20, 210, 100, 20, hwnd, NULL, NULL, NULL);
        hGridSpacing = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "10", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                      150, 210, 100, 25, hwnd, (HMENU)IDC_GRID_SPACING, NULL, NULL);
        CreateWindow("STATIC",
                     "Spacing between grid nodes for TPS interpolation.",
                     WS_VISIBLE | WS_CHILD, 260, 210, 500, 40, hwnd, NULL, NULL, NULL);

        // Max TPS Points Field
        CreateWindowEx(WS_EX_CLIENTEDGE, "STATIC", "Max TPS Points:", WS_VISIBLE | WS_CHILD,
                       20, 260, 120, 20, hwnd, NULL, NULL, NULL);
        hMaxTPSPoints = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "0", WS_CHILD | WS_VISIBLE | WS_BORDER,
                                       150, 260, 100, 25, hwnd, (HMENU)IDC_MAX_TPS_POINTS, NULL, NULL);
        CreateWindow("STATIC",
                     "Maximum number of points for TPS computation (0 = use all).",
                     WS_VISIBLE | WS_CHILD, 260, 260, 500, 40, hwnd, NULL, NULL, NULL);

        // Run Button
        CreateWindow("BUTTON", "Run xyz2dxf", WS_VISIBLE | WS_CHILD | BS_PUSHBUTTON,
                     20, 320, 200, 30, hwnd, (HMENU)IDC_RUN_BUTTON, NULL, NULL);

        break;
    }
    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_BROWSE_BUTTON)
        {
            std::string filePath = openFileDialog(hwnd);
            if (!filePath.empty())
            {
                SetWindowText(hInputFile, filePath.c_str());
            }
        }
        else if (LOWORD(wParam) == IDC_RUN_BUTTON)
        {
            // Collect inputs and build the command
            char inputFile[260], minDist[64], precision[64], pdMode[64], gridSpacing[64], maxTPSPoints[64];
            GetWindowText(hInputFile, inputFile, sizeof(inputFile));
            GetWindowText(hMinDist, minDist, sizeof(minDist));
            GetWindowText(hPrecision, precision, sizeof(precision));
            GetWindowText(hPDMODE, pdMode, sizeof(pdMode));
            GetWindowText(hGridSpacing, gridSpacing, sizeof(gridSpacing));
            GetWindowText(hMaxTPSPoints, maxTPSPoints, sizeof(maxTPSPoints));

            std::stringstream cmd;
            cmd << "xyz2dxf " << "\"" << inputFile << "\" " << minDist << " " << precision << " " << pdMode;
            if (strlen(gridSpacing) > 0)
                cmd << " " << gridSpacing;
            if (strlen(maxTPSPoints) > 0)
                cmd << " " << maxTPSPoints;

            std::string command = cmd.str();

            // Execute the command
            STARTUPINFO si = {};
            PROCESS_INFORMATION pi = {};
            si.cb = sizeof(si);

            if (!CreateProcess(NULL, const_cast<char *>(command.c_str()), NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi))
            {
                MessageBox(hwnd, "Failed to run xyz2dxf.exe. Ensure it's in the PATH.", "Error", MB_ICONERROR | MB_OK);
            }
            else
            {
                CloseHandle(pi.hProcess);
                CloseHandle(pi.hThread);
            }
        }
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

/**
 * WinMain:
 * --------
 * Entry point for the Windows application. Registers the window class, creates the main window,
 * and runs the message loop.
 */
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    const char CLASS_NAME[] = "xyz2dxfWindowClass";

    // Define and register the window class
    WNDCLASS wc = {};
    wc.lpfnWndProc = WindowProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;

    RegisterClass(&wc);

    // Create the main application window
    HWND hwnd = CreateWindowEx(0, CLASS_NAME, "XYZ to DXF Converter", WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU,
                               CW_USEDEFAULT, CW_USEDEFAULT, 800, 400, NULL, NULL, hInstance, NULL);

    if (!hwnd)
        return 0;

    ShowWindow(hwnd, nCmdShow);

    // Message loop
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return 0;
}
