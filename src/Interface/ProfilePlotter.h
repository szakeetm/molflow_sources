
#pragma once

#include "GLApp/GLWindow.h"
#include "GLApp/GLChart/GLChartConst.h"
#include <vector>
#include <map>

class GLChart;
class GLLabel;
class GLCombo;
class GLButton;
class GLFormula;
class GLDataView;
class GLToggle;
class GLTextField;
class Worker;
class InterfaceGeometry;

class ProfilePlotter : public GLWindow {

public:

    // Construction
    ProfilePlotter(Worker* work);

    // Component method
    void Display(Worker* w);
    void Refresh();
    void Update(float appTime, bool force = false);
    void Reset();
    void RefreshPlottedColors();

    // Implementation
    void ProcessMessage(GLComponent* src, int message) override;
    void SetBounds(int x, int y, int w, int h) override;
    int addView(int facet);
    std::vector<int> GetViews();
    void SetViews(const std::vector<int>& updatedViews);
    bool IsLogScaled();
    void SetLogScaled(bool logScale);
    void SetWorker(Worker* w);
    std::map<int, GLColor> GetIDColorPairs() const;

private:
    int remView(int facet);
    void refreshViews();
    void applyFacetHighlighting();

    Worker* worker;
    GLButton* dismissButton;
    GLChart* chart;
    GLCombo* profCombo;
    GLTextField* selFacInput;
    GLLabel* normLabel;
    GLLabel* warningLabel;
    GLCombo* displayModeCombo;
    //GLToggle    *showAllMoments;

    GLButton* selButton;
    GLButton* addButton;
    GLButton* removeButton;
    GLButton* removeAllButton;
    GLTextField* formulaText;
    GLButton* plotExpressionBtn;
    GLToggle* logYToggle;
    GLToggle* correctForGas;

    GLToggle* colorToggle;
    GLLabel* fixedLineWidthText;
    GLButton* fixedLineWidthButton;
    GLTextField* fixedLineWidthField;
    GLToggle* useProfColToggle;
    GLButton* selectPlottedButton;

    GLDataView* views[MAX_VIEWS];
    int          nbView;
    float        lastUpdate;
    std::map<int, GLColor> plottedFacets;
};