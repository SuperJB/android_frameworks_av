LOCAL_PATH:= $(call my-dir)
include $(CLEAR_VARS)

LOCAL_SRC_FILES:= \
	Camera.cpp \
	CameraParameters.cpp \
	ICamera.cpp \
	ICameraClient.cpp \
	ICameraService.cpp \
	ICameraRecordingProxy.cpp \
	ICameraRecordingProxyListener.cpp

LOCAL_SHARED_LIBRARIES := \
	libcutils \
	libutils \
	libbinder \
	libhardware \
	libui \
	libgui

ifeq ($(TARGET_BOARD_PLATFORM),exDroid)
	LOCAL_CFLAGS += -DALLWINNER
endif

ifeq ($(BOARD_USES_QCOM_HARDWARE),true)
	LOCAL_CFLAGS += -DQCOM_HARDWARE
endif

ifeq ($(BOARD_CAMERA_HAVE_ISO),true)
	LOCAL_CFLAGS += -DHAVE_ISO
endif

LOCAL_MODULE:= libcamera_client

include $(BUILD_SHARED_LIBRARY)
