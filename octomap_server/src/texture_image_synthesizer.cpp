#include <ros/ros.h>
#include <octomap_msgs/conversions.h>
#include <octomap/octomap.h>
#include <fstream>
#include <cv_bridge/cv_bridge.h>

#include <octomap_msgs/TextureViewSynthesis.h>
using octomap_msgs::TextureViewSynthesis;

#include <geometry_msgs/Pose.h>
#include <sensor_msgs/Image.h>

using namespace std;
using namespace octomap;

class TextureImageSynthesizer{
public:
  TextureImageSynthesizer(){
   
    image_pub = n.advertise<sensor_msgs::Image>("image", 1, true);
    depth_pub = n.advertise<sensor_msgs::Image>("depth", 1, true);
    pose_sub = n.subscribe("pose", 1, &TextureImageSynthesizer::poseCallback, this);
  }
protected:
  ros::Subscriber pose_sub;
  ros::Publisher image_pub, depth_pub;
  ros::NodeHandle n;
  
  void poseCallback(const geometry_msgs::Pose& msg) 
  {
    //ros::Time start = ros::Time::now();
    std::vector<geometry_msgs::Pose> poses;
    for(unsigned i = 0; i<10; ++i)
      poses.push_back(msg);
    TextureViewSynthesis::Request req;
    req.poses = poses;
    req.w = 160;
    req.h = 120;
    req.fx = 100;
    req.fy = 100;
    req.cx = 80;
    req.cy = 60;
    std::string servname = "/octomap_server/synthesize_views";
    ROS_INFO("Sending view synthesis request to %s...", n.resolveName(servname).c_str());
    TextureViewSynthesis::Response resp;
    while(n.ok() && !ros::service::call(servname, req, resp))
    {
      ROS_WARN("Request to %s failed; trying again...", n.resolveName(servname).c_str());
      usleep(1000000);
    }

    image_pub.publish(resp.images.at(0));
    depth_pub.publish(resp.depths.at(0));
    //ros::Time end = ros::Time::now();
    //double elapsed = (end-start).toSec();
    //ROS_INFO("Time to synthesize one view: %f", elapsed);
  }
};

int main(int argc, char** argv){
  ros::init(argc, argv, "texture_image_synthesizer");
  try{
    TextureImageSynthesizer tis;
    ros::spin();
  }catch(std::runtime_error& e){
    ROS_ERROR("texture_image_synthesizer exception: %s", e.what());
    exit(2);
  }

  exit(0);
}


