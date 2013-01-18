LOCAL_PATH:= $(call my-dir)

#
# libmediaplayerservice
#

include $(CLEAR_VARS)

ifeq ($(TARGET_BOARD_PLATFORM),exDroid)
	LOCAL_CFLAGS += -DALLWINNER
endif

LOCAL_SRC_FILES:=               \
    ActivityManager.cpp         \
    Crypto.cpp                  \
    HDCP.cpp                    \
    MediaPlayerFactory.cpp      \
    MediaPlayerService.cpp      \
    MediaRecorderClient.cpp     \
    MetadataRetrieverClient.cpp \
    MidiFile.cpp                \
    MidiMetadataRetriever.cpp   \
    RemoteDisplay.cpp           \
    StagefrightPlayer.cpp       \
    StagefrightRecorder.cpp

LOCAL_SHARED_LIBRARIES :=     		\
	libcutils             			\
	libutils              			\
	libbinder             			\
	libvorbisidec         			\
	libsonivox            			\
	libmedia              			\
	libmedia_native       			\
	libcamera_client      			\
	libstagefright        			\
	libstagefright_omx    			\
	libstagefright_foundation       \
	libgui                          \
	libdl                           \
    libstagefright_wfd

ifeq ($(TARGET_BOARD_PLATFORM),exDroid)
LOCAL_SRC_FILES+= \
    CedarPlayer.cpp       		\
    CedarAPlayerWrapper.cpp		\
    SimpleMediaFormatProbe.cpp	\
    MovAvInfoDetect.cpp

LOCAL_SHARED_LIBRARIES += \
	libCedarX           			\
	libCedarA
endif

LOCAL_STATIC_LIBRARIES := \
        libstagefright_nuplayer                 \
        libstagefright_rtsp                     \

LOCAL_C_INCLUDES :=                                                 \
    $(call include-path-for, graphics corecg)                       \
    $(TOP)/frameworks/av/media/libstagefright/include               \
    $(TOP)/frameworks/av/media/libstagefright/rtsp                  \
    $(TOP)/frameworks/av/media/libstagefright/wifi-display          \
    $(TOP)/frameworks/native/include/media/openmax                  \
    $(TOP)/external/tremolo/Tremolo                                 \

ifeq ($(TARGET_BOARD_PLATFORM),exDroid)
LOCAL_C_INCLUDES +=  \
        $(TOP)/external/cedarx/CedarXAndroid/IceCreamSandwich \
        $(TOP)/external/cedarx/CedarX/include/include_audio \
        $(TOP)/external/cedarx/CedarX/include/include_cedarv \
        $(TOP)/external/cedarx/CedarX/include \
        $(TOP)/external/cedarx/CedarA \
        $(TOP)/external/cedarx/CedarA/include 
endif

ifeq ($(TARGET_BOARD_PLATFORM),exDroid)
	LOCAL_CFLAGS +=-DCEDARX_ANDROID_VERSION=7
endif

ifeq ($(BOARD_USES_QCOM_HARDWARE),true)
    LOCAL_STATIC_LIBRARIES += libmedia_helper
	LOCAL_C_INCLUDES += $(TOP)/frameworks/av/include/media
endif

LOCAL_MODULE:= libmediaplayerservice

include $(BUILD_SHARED_LIBRARY)

include $(call all-makefiles-under,$(LOCAL_PATH))
