
#pragma once

#include "GlobalSettings_shared.h"

class GlobalSettings : public GlobalSettingsBase {

public:

  // Construction
  GlobalSettings(Worker *w);
  void UpdateOutgassing();

  // Implementation
  void ProcessMessage(GLComponent *src,int message);
  void Update();
  void ResizeProcessPanel(int windowWidth, int windowHeight);
  void SetBounds(int x, int y, int width, int height) override;

private:

  GLTextField *outgassingGasRateText;
  GLTextField *outgassingMoleculeRateText;
  GLLabel     *desorbedMoleculesLabel;
  GLTextField *desorbedMoleculesText;
  GLTextField *gasMassText;
  GLToggle    *enableDecay;
  GLTextField *halfLifeText;
  GLButton* recalcButton;
  GLToggle* useOldXMLFormat;

};
