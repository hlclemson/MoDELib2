/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalActor_h_
#define model_DefectiveCrystalActor_h_

#include <QGroupBox>
#include <QGridLayout>
#include <QCheckBox>
#include <QLabel>
#include <QTabWidget>
#include <QTextEdit>
#include <QLineEdit>
#include <QPushButton>
#include <QVTKOpenGLStereoWidget.h>

#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <DislocationDynamicsBase.h>
//#include <ConfigurationFields.h>
#include <DefectiveCrystal.h>
#include <DDauxIO.h>
#include <DDFieldWidget.h>
#include <SimplicialMesh.h>
#include <NetworkNodeActor.h>
#include <NetworkLinkActor.h>
#include <NetworkLoopActor.h>
#include <InclusionActor.h>
#include <GlidePlaneActor.h>
#include <QuadratureActor.h>
#include <ChartActor.h>
#include <QCheckBox>

// added std library
#include <atomic>

namespace model
{
    struct DefectiveCrystalActor : public QWidget
                        ,private DDconfigIO<3>
                        ,private DDauxIO<3>
    {
        
        Q_OBJECT
        
        
    public:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        QVTKOpenGLStereoWidget* const qvtkGLwidget;
//        const DDtraitsIO& traitsIO;
        DislocationDynamicsBase<3>& ddBase;
        DefectiveCrystal<3> defectiveCrystal;
        NetworkNodeActor* nodes;
        NetworkLinkActor* segments;
        NetworkLoopActor* loops;
        InclusionActor* inclusions;
        GlidePlaneActor* glidePlanes;
        QuadratureActor* quadrature;
        ChartActor* chartActor;
        DDFieldWidget* ddField;
    private:
        std::atomic<bool> m_stopRequested{false};

        private slots:
        bool updateConfiguration();
        bool nextConfiguration();
        bool prevConfiguration();
        void playConfigurations();
        void stopConfigurations();

        public:
        
        QGridLayout* mainLayout;
        QLineEdit* frameIDedit;
        QPushButton* plusFrameButton;
        QPushButton* playFrameButton;
        QPushButton* stopFrameButton;
        QPushButton* minusFrameButton;
        QLineEdit* frameIncrementEdit;
        QCheckBox* saveImage;
        QTabWidget* tabWidget;

        DefectiveCrystalActor(vtkGenericOpenGLRenderWindow* const,
                    vtkRenderer* const ren,
                    QVTKOpenGLStereoWidget* const qvtkGLwidget_in,
                    DislocationDynamicsBase<3>& ddBase_in);
        bool updateConfiguration(const size_t& frameID);
        void modify();
        const DDconfigIO<3>& configIO() const;
        DDconfigIO<3>& configIO();
        const DDauxIO<3>& auxIO() const;
        DDauxIO<3>& auxIO();

    };
    
}
#endif
