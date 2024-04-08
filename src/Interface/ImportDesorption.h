
#include "GLApp/GLWindow.h"
#include <vector>
class GLButton;
class GLTextField;
class GLLabel;
class GLToggle;
class GLTitledPanel;
class GLCombo;

class InterfaceGeometry;
class Worker;

#ifndef _IMPORTDESORPTIONH_
#define _IMPORTDESORPTIONH_

class ImportDesorption : public GLWindow {

public:

  // Construction
  ImportDesorption();

  // Component methods
  void SetGeometry(InterfaceGeometry *s,Worker *w);

  // Implementation
  void ProcessMessage(GLComponent *src,int message) override;

private:

	void    LoadConvFile(const char* fileName);
	void   EnableDisableComponents();

  InterfaceGeometry      *interfGeom;
  Worker	    *work;

  GLTitledPanel *filePanel;
  GLTextField   *synFileName;
  GLButton      *loadSynButton,*reloadButton;
  GLButton      *useCurrentButton;
  GLLabel       *analysisResultLabel;

  GLTitledPanel *importPanel;
  GLCombo       *sourceSelectCombo;
  GLTextField   *timeField;

  GLToggle      *r1,*r2,*r3;
  GLTextField   *eta0Text,*alphaText,*cutoffText;
  GLTextField   *convFileName;
  GLButton      *loadConvButton;
  GLButton      *convInfoButton;
  GLLabel       *convAnalysisLabel;

  GLButton      *setButton;
  GLButton      *cancelButton;

  size_t doseSource;
  size_t mode;
  std::vector<std::pair<double,double>> convDistr;
};

#endif /* _IMPORTDESORPTIONH_ */

