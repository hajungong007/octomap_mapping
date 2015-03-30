/*
 * Copyright (c) 2010-2013, A. Hornung, University of Freiburg
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <octomap_server/TextureOctomapServer.h>

using namespace octomap;
using octomap_msgs::Octomap;

namespace octomap_server{

TextureOctomapServer::TextureOctomapServer(ros::NodeHandle private_nh_)
: m_nh(),
  m_pointCloudSub(NULL),
  m_tfPointCloudSub(NULL),
  m_octree(NULL),
  m_maxRange(-1.0),
  m_stereoErrorCoeff(10.0),
  m_worldFrameId("/map"),
  m_useHeightMap(true),
  m_colorFactor(0.8),
  m_latchedTopics(true),
  m_publishFreeSpace(false),
  m_res(0.05),
  m_treeDepth(0),
  m_maxTreeDepth(0),
  m_occThres(0.5),
  m_probHit(0.7), m_probMiss(0.4),
  m_thresMin(0.12), m_thresMax(0.97),
  m_pointcloudMinZ(-std::numeric_limits<double>::max()),
  m_pointcloudMaxZ(std::numeric_limits<double>::max()),
  m_occupancyMinZ(-std::numeric_limits<double>::max()),
  m_occupancyMaxZ(std::numeric_limits<double>::max()),
  m_minSizeX(0.0), m_minSizeY(0.0),
  m_filterSpeckles(false),
  m_compressMap(true),
  m_incrementalUpdate(false)
{
  ros::NodeHandle private_nh(private_nh_);
  private_nh.param("frame_id", m_worldFrameId, m_worldFrameId);
  private_nh.param("height_map", m_useHeightMap, m_useHeightMap);
  private_nh.param("color_factor", m_colorFactor, m_colorFactor);

  private_nh.param("pointcloud_min_z", m_pointcloudMinZ,m_pointcloudMinZ);
  private_nh.param("pointcloud_max_z", m_pointcloudMaxZ,m_pointcloudMaxZ);
  private_nh.param("occupancy_min_z", m_occupancyMinZ,m_occupancyMinZ);
  private_nh.param("occupancy_max_z", m_occupancyMaxZ,m_occupancyMaxZ);
  private_nh.param("min_x_size", m_minSizeX,m_minSizeX);
  private_nh.param("min_y_size", m_minSizeY,m_minSizeY);

  private_nh.param("filter_speckles", m_filterSpeckles, m_filterSpeckles);

  private_nh.param("sensor_model/max_range", m_maxRange, m_maxRange);
  private_nh.param("sensor_model/stereo_error_coeff", m_stereoErrorCoeff, m_stereoErrorCoeff);

  private_nh.param("resolution", m_res, m_res);
  private_nh.param("occupancy_threshold", m_occThres, m_occThres);
  private_nh.param("sensor_model/hit", m_probHit, m_probHit);
  private_nh.param("sensor_model/miss", m_probMiss, m_probMiss);
  private_nh.param("sensor_model/min", m_thresMin, m_thresMin);
  private_nh.param("sensor_model/max", m_thresMax, m_thresMax);
  private_nh.param("compress_map", m_compressMap, m_compressMap);
  private_nh.param("incremental_2D_projection", m_incrementalUpdate, m_incrementalUpdate);


  // initialize octomap object & params
  ROS_INFO("OctoMap Server: Using stereo error model.");
  m_octree = new TextureOcTree(m_res, m_stereoErrorCoeff, m_maxRange);
  m_octree->setOccupancyThres(m_occThres);
  m_octree->setProbHit(m_probHit);
  m_octree->setProbMiss(m_probMiss);
  m_octree->setClampingThresMin(m_thresMin);
  m_octree->setClampingThresMax(m_thresMax);
  m_treeDepth = m_octree->getTreeDepth();
  m_maxTreeDepth = m_treeDepth;
  m_gridmap.info.resolution = m_res;

  double r, g, b, a;
  private_nh.param("color/r", r, 0.0);
  private_nh.param("color/g", g, 0.0);
  private_nh.param("color/b", b, 1.0);
  private_nh.param("color/a", a, 1.0);
  m_color.r = r;
  m_color.g = g;
  m_color.b = b;
  m_color.a = a;

  private_nh.param("color_free/r", r, 0.0);
  private_nh.param("color_free/g", g, 1.0);
  private_nh.param("color_free/b", b, 0.0);
  private_nh.param("color_free/a", a, 1.0);
  m_colorFree.r = r;
  m_colorFree.g = g;
  m_colorFree.b = b;
  m_colorFree.a = a;

  private_nh.param("publish_free_space", m_publishFreeSpace, m_publishFreeSpace);

  private_nh.param("latch", m_latchedTopics, m_latchedTopics);
  if (m_latchedTopics){
    ROS_INFO("Publishing latched (single publish will take longer, all topics are prepared)");
  } else
    ROS_INFO("Publishing non-latched (topics are only prepared as needed, will only be re-published on map change");

  m_markerPub = m_nh.advertise<visualization_msgs::MarkerArray>("occupied_cells_vis_array", 1, m_latchedTopics);
  m_binaryMapPub = m_nh.advertise<Octomap>("octomap_binary", 1, m_latchedTopics);
  m_fullMapPub = m_nh.advertise<Octomap>("octomap_full", 1, m_latchedTopics);
  m_pointCloudPub = m_nh.advertise<sensor_msgs::PointCloud2>("octomap_point_cloud_centers", 1, m_latchedTopics);
  m_mapPub = m_nh.advertise<nav_msgs::OccupancyGrid>("projected_map", 5, m_latchedTopics);	
  m_fmarkerPub = m_nh.advertise<visualization_msgs::MarkerArray>("free_cells_vis_array", 1, m_latchedTopics);

  m_pointCloudSub = new message_filters::Subscriber<sensor_msgs::PointCloud2> (m_nh, "cloud_in", 5);
  ROS_INFO("Setting up message filter with world frame %s", m_worldFrameId.c_str());
  m_tfPointCloudSub = new tf::MessageFilter<sensor_msgs::PointCloud2> (*m_pointCloudSub, m_tfListener, m_worldFrameId, 5);
  m_tfPointCloudSub->registerCallback(boost::bind(&TextureOctomapServer::insertCloudCallback, this, _1));

  m_octomapBinaryService = m_nh.advertiseService("octomap_binary", &TextureOctomapServer::octomapBinarySrv, this);
  m_octomapFullService = m_nh.advertiseService("octomap_full", &TextureOctomapServer::octomapFullSrv, this);
  m_clearBBXService = private_nh.advertiseService("clear_bbx", &TextureOctomapServer::clearBBXSrv, this);
  m_resetService = private_nh.advertiseService("reset", &TextureOctomapServer::resetSrv, this);
  m_viewSynthesisService = private_nh.advertiseService("synthesize_views", &TextureOctomapServer::synthesizeViewsSrv, this);

  dynamic_reconfigure::Server<OctomapServerConfig>::CallbackType f;

  f = boost::bind(&TextureOctomapServer::reconfigureCallback, this, _1, _2);
  m_reconfigureServer.setCallback(f);
  
  // Set up change tracking
  std::string changeSetTopic = "changes";
  private_nh.param("topic_changes", changeSetTopic, changeSetTopic);
  private_nh.param("track_changes", track_changes, false);

  if (track_changes) {
    ROS_INFO("Octmap server set to track changes");
    pubChangeSet = private_nh.advertise<sensor_msgs::PointCloud2>(changeSetTopic, 1);
    m_octree->enableChangeDetection(true);
  }
}

TextureOctomapServer::~TextureOctomapServer(){
  if (m_tfPointCloudSub){
    delete m_tfPointCloudSub;
    m_tfPointCloudSub = NULL;
  }

  if (m_pointCloudSub){
    delete m_pointCloudSub;
    m_pointCloudSub = NULL;
  }


  if (m_octree){
    delete m_octree;
    m_octree = NULL;
  }

}

bool TextureOctomapServer::openFile(const std::string& filename){
  if (filename.length() <= 3)
    return false;

  std::string suffix = filename.substr(filename.length()-3, 3);
  if (suffix== ".bt"){
    if (!m_octree->readBinary(filename)){
      return false;
    }
  } else if (suffix == ".ot"){
    AbstractOcTree* tree = AbstractOcTree::read(filename);
    if (!tree){
      return false;
    }
    if (m_octree){
      delete m_octree;
      m_octree = NULL;
    }
    m_octree = dynamic_cast<TextureOcTree*>(tree);
    if (!m_octree){
      ROS_ERROR("Could not read TextureOcTree in file, currently there are no other types supported in .ot");
      return false;
    }

  } else{
    return false;
  }

  ROS_INFO("Octomap file %s loaded (%zu nodes).", filename.c_str(),m_octree->size());

  m_treeDepth = m_octree->getTreeDepth();
  m_maxTreeDepth = m_treeDepth;
  m_res = m_octree->getResolution();
  m_gridmap.info.resolution = m_res;
  double minX, minY, minZ;
  double maxX, maxY, maxZ;
  m_octree->getMetricMin(minX, minY, minZ);
  m_octree->getMetricMax(maxX, maxY, maxZ);

  m_updateBBXMin[0] = m_octree->coordToKey(minX);
  m_updateBBXMin[1] = m_octree->coordToKey(minY);
  m_updateBBXMin[2] = m_octree->coordToKey(minZ);
  
  m_updateBBXMax[0] = m_octree->coordToKey(maxX);
  m_updateBBXMax[1] = m_octree->coordToKey(maxY);
  m_updateBBXMax[2] = m_octree->coordToKey(maxZ);
  
  publishAll();

  return true;

}

void TextureOctomapServer::insertCloudCallback(const sensor_msgs::PointCloud2::ConstPtr& cloud){
  ros::WallTime startTime = ros::WallTime::now();

  PCLPointCloudRGB pc_rgb; // input cloud for filtering
  pcl::fromROSMsg(*cloud, pc_rgb);
  PCLPointCloud pc;
  pcl::PointCloudXYZRGBtoXYZI(pc_rgb, pc);

  tf::StampedTransform sensorToWorldTf;
  try {
    m_tfListener.lookupTransform(m_worldFrameId, cloud->header.frame_id, cloud->header.stamp, sensorToWorldTf);
  } catch(tf::TransformException& ex){
    ROS_ERROR_STREAM( "Transform error of sensor data: " << ex.what() << ", quitting callback");
    return;
  }

  Eigen::Matrix4f sensorToWorld;
  pcl_ros::transformAsMatrix(sensorToWorldTf, sensorToWorld);

  // Get pose of sensor as vector for orientation of +x for sensor frame
  // First get rotation from world to sensor frame
  Eigen::Matrix3f worldToSensorRot = (sensorToWorld.inverse()).block<3,3>(0,0);
  Eigen::Vector3f sensorOrientation = worldToSensorRot * Eigen::Vector3f(1.0,0.0,0.0);
  octomath::Vector3 orientation(sensorOrientation(0), sensorOrientation(1), sensorOrientation(2)); 


  // set up filter for height range, also removes NANs:
  pcl::PassThrough<pcl::PointXYZI> pass;
  pass.setFilterFieldName("z");
  pass.setFilterLimits(m_pointcloudMinZ, m_pointcloudMaxZ);

  // directly transform to map frame:
  pcl::transformPointCloud(pc, pc, sensorToWorld);

  // just filter height range:
  pass.setInputCloud(pc.makeShared());
  pass.filter(pc);

  insertScan(sensorToWorldTf.getOrigin(), orientation, pc);

  double total_elapsed = (ros::WallTime::now() - startTime).toSec();
  ROS_DEBUG("Pointcloud insertion in TextureOctomapServer done (%zu pts, %f sec)", pc.size(), total_elapsed);

  publishAll(cloud->header.stamp);
}

void TextureOctomapServer::insertScan(const tf::Point& sensorOriginTf, 
                                      const octomath::Vector3& sensorOrientation,
                                      const PCLPointCloud& scan)
{
  point3d sensorOrigin = pointTfToOctomap(sensorOriginTf);

  if (!m_octree->coordToKeyChecked(sensorOrigin, m_updateBBXMin)
    || !m_octree->coordToKeyChecked(sensorOrigin, m_updateBBXMax))
  {
    ROS_ERROR_STREAM("Could not generate Key for origin "<<sensorOrigin);
  }

  // Insert points: free on ray, occupied on endpoint:
  KeySet free_cells, occupied_cells;
  for (PCLPointCloud::const_iterator it = scan.begin(); it != scan.end(); ++it){
    point3d point(it->x, it->y, it->z);
    // maxrange check
    if ((m_maxRange < 0.0) || ((point - sensorOrigin).norm() <= m_maxRange) ) {

      // free cells
      if (m_octree->computeRayKeys(sensorOrigin, point, m_keyRay)){
        free_cells.insert(m_keyRay.begin(), m_keyRay.end());
      }
      // occupied endpoint
      OcTreeKey key;
      if (m_octree->coordToKeyChecked(point, key)){
        occupied_cells.insert(key);

        updateMinKey(key, m_updateBBXMin);
        updateMaxKey(key, m_updateBBXMax);
      }
    } else {// ray longer than maxrange:;
      point3d new_end = sensorOrigin + (point - sensorOrigin).normalized() * m_maxRange;
      if (m_octree->computeRayKeys(sensorOrigin, new_end, m_keyRay)){
        free_cells.insert(m_keyRay.begin(), m_keyRay.end());

        octomap::OcTreeKey endKey;
        if (m_octree->coordToKeyChecked(new_end, endKey)){
          updateMinKey(endKey, m_updateBBXMin);
          updateMaxKey(endKey, m_updateBBXMax);
        } else{
          ROS_ERROR_STREAM("Could not generate Key for endpoint "<<new_end);
        }
      }
    }
  }

  // insert data into tree  -----------------------
  for (KeySet::iterator it = free_cells.begin(); it != free_cells.end(); ++it) {
    float d = fabs(((m_octree->keyToCoord(*it)) - sensorOrigin).dot(sensorOrientation));
    m_octree->updateNode(*it, false, d);
  }
  for (KeySet::iterator it = occupied_cells.begin(); it != occupied_cells.end(); ++it) {
    float d = fabs(((m_octree->keyToCoord(*it)) - sensorOrigin).dot(sensorOrientation));
    m_octree->updateNode(*it, true, d);
  }

  // Now update intensities from all points not filtered
  for (PCLPointCloud::const_iterator it = scan.begin(); it != scan.end(); ++it)
  {
    point3d point(it->x, it->y, it->z);
    unsigned char intensity = floor(it->intensity);
    m_octree->insertTexturePoint(point,intensity,sensorOrigin);
  }

  if (track_changes) {
    trackChanges();
  }

  // TODO: eval lazy+updateInner vs. proper insertion
  // non-lazy by default (updateInnerOccupancy() too slow for large maps)
  //m_octree->updateInnerOccupancy();
  octomap::point3d minPt, maxPt;
  ROS_DEBUG_STREAM("Bounding box keys (before): " << m_updateBBXMin[0] << " " <<m_updateBBXMin[1] << " " << m_updateBBXMin[2] << " / " <<m_updateBBXMax[0] << " "<<m_updateBBXMax[1] << " "<< m_updateBBXMax[2]);


  // TODO: we could also limit the bbx to be within the map bounds here (see publishing check)
  minPt = m_octree->keyToCoord(m_updateBBXMin);
  maxPt = m_octree->keyToCoord(m_updateBBXMax);
  ROS_DEBUG_STREAM("Updated area bounding box: "<< minPt << " - "<<maxPt);
  ROS_DEBUG_STREAM("Bounding box keys (after): " << m_updateBBXMin[0] << " " <<m_updateBBXMin[1] << " " << m_updateBBXMin[2] << " / " <<m_updateBBXMax[0] << " "<<m_updateBBXMax[1] << " "<< m_updateBBXMax[2]);

  if (m_compressMap)
    m_octree->prune();
}



void TextureOctomapServer::publishAll(const ros::Time& rostime){
  ros::WallTime startTime = ros::WallTime::now();
  size_t octomapSize = m_octree->size();
  // TODO: estimate num occ. voxels for size of arrays (reserve)
  if (octomapSize <= 1){
    ROS_WARN("Nothing to publish, octree is empty");
    return;
  }

  bool publishFreeMarkerArray = m_publishFreeSpace && (m_latchedTopics || m_fmarkerPub.getNumSubscribers() > 0);
  bool publishMarkerArray = (m_latchedTopics || m_markerPub.getNumSubscribers() > 0);
  bool publishPointCloud = (m_latchedTopics || m_pointCloudPub.getNumSubscribers() > 0);
  bool publishBinaryMap = (m_latchedTopics || m_binaryMapPub.getNumSubscribers() > 0);
  bool publishFullMap = (m_latchedTopics || m_fullMapPub.getNumSubscribers() > 0);
  m_publish2DMap = (m_latchedTopics || m_mapPub.getNumSubscribers() > 0);

  // init markers for free space:
  visualization_msgs::MarkerArray freeNodesVis;
  // each array stores all cubes of a different size, one for each depth level:
  freeNodesVis.markers.resize(m_treeDepth+1);

  geometry_msgs::Pose pose;
  pose.orientation = tf::createQuaternionMsgFromYaw(0.0);

  // init markers:
  visualization_msgs::MarkerArray occupiedNodesVis;
  // each array stores all cubes of a different size, one for each depth level:
  occupiedNodesVis.markers.resize(m_treeDepth+1);

  // init pointcloud:
  pcl::PointCloud<pcl::PointXYZ> pclCloud;

  // call pre-traversal hook:
  handlePreNodeTraversal(rostime);

  // now, traverse all leafs in the tree:
  for (OcTreeT::iterator it = m_octree->begin(m_maxTreeDepth),
      end = m_octree->end(); it != end; ++it)
  {
    bool inUpdateBBX = isInUpdateBBX(it);

    // call general hook:
    handleNode(it);
    if (inUpdateBBX)
      handleNodeInBBX(it);

    if (m_octree->isNodeOccupied(*it)){
      double z = it.getZ();
      if (z > m_occupancyMinZ && z < m_occupancyMaxZ)
      {
        double size = it.getSize();
        double x = it.getX();
        double y = it.getY();

        // Ignore speckles in the map:
        if (m_filterSpeckles && (it.getDepth() == m_treeDepth +1) && isSpeckleNode(it.getKey())){
          ROS_DEBUG("Ignoring single speckle at (%f,%f,%f)", x, y, z);
          continue;
        } // else: current octree node is no speckle, send it out

        handleOccupiedNode(it);
        if (inUpdateBBX)
          handleOccupiedNodeInBBX(it);


        //create marker:
        if (publishMarkerArray){
          unsigned idx = it.getDepth();
          assert(idx < occupiedNodesVis.markers.size());

          geometry_msgs::Point cubeCenter;
          cubeCenter.x = x;
          cubeCenter.y = y;
          cubeCenter.z = z;

          occupiedNodesVis.markers[idx].points.push_back(cubeCenter);
          if (m_useHeightMap){
            double minX, minY, minZ, maxX, maxY, maxZ;
            m_octree->getMetricMin(minX, minY, minZ);
            m_octree->getMetricMax(maxX, maxY, maxZ);

            double h = (1.0 - std::min(std::max((cubeCenter.z-minZ)/ (maxZ - minZ), 0.0), 1.0)) *m_colorFactor;
            occupiedNodesVis.markers[idx].colors.push_back(heightMapColor(h));
          }
        }

        // insert into pointcloud:
        if (publishPointCloud)
          pclCloud.push_back(pcl::PointXYZ(x, y, z));

      }
    } else{ // node not occupied => mark as free in 2D map if unknown so far
      double z = it.getZ();
      if (z > m_occupancyMinZ && z < m_occupancyMaxZ)
      {
        handleFreeNode(it);
        if (inUpdateBBX)
          handleFreeNodeInBBX(it);

        if (m_publishFreeSpace){
          double x = it.getX();
          double y = it.getY();

          //create marker for free space:
          if (publishFreeMarkerArray){
            unsigned idx = it.getDepth();
            assert(idx < freeNodesVis.markers.size());

            geometry_msgs::Point cubeCenter;
            cubeCenter.x = x;
            cubeCenter.y = y;
            cubeCenter.z = z;

            freeNodesVis.markers[idx].points.push_back(cubeCenter);
          }
        }

      }
    }
  }

  // call post-traversal hook:
  handlePostNodeTraversal(rostime);

  // finish MarkerArray:
  if (publishMarkerArray){
    for (unsigned i= 0; i < occupiedNodesVis.markers.size(); ++i){
      double size = m_octree->getNodeSize(i);

      occupiedNodesVis.markers[i].header.frame_id = m_worldFrameId;
      occupiedNodesVis.markers[i].header.stamp = rostime;
      occupiedNodesVis.markers[i].ns = "map";
      occupiedNodesVis.markers[i].id = i;
      occupiedNodesVis.markers[i].type = visualization_msgs::Marker::CUBE_LIST;
      occupiedNodesVis.markers[i].scale.x = size;
      occupiedNodesVis.markers[i].scale.y = size;
      occupiedNodesVis.markers[i].scale.z = size;
      occupiedNodesVis.markers[i].color = m_color;


      if (occupiedNodesVis.markers[i].points.size() > 0)
        occupiedNodesVis.markers[i].action = visualization_msgs::Marker::ADD;
      else
        occupiedNodesVis.markers[i].action = visualization_msgs::Marker::DELETE;
    }

    m_markerPub.publish(occupiedNodesVis);
  }


  // finish FreeMarkerArray:
  if (publishFreeMarkerArray){
    for (unsigned i= 0; i < freeNodesVis.markers.size(); ++i){
      double size = m_octree->getNodeSize(i);

      freeNodesVis.markers[i].header.frame_id = m_worldFrameId;
      freeNodesVis.markers[i].header.stamp = rostime;
      freeNodesVis.markers[i].ns = "map";
      freeNodesVis.markers[i].id = i;
      freeNodesVis.markers[i].type = visualization_msgs::Marker::CUBE_LIST;
      freeNodesVis.markers[i].scale.x = size;
      freeNodesVis.markers[i].scale.y = size;
      freeNodesVis.markers[i].scale.z = size;
      freeNodesVis.markers[i].color = m_colorFree;


      if (freeNodesVis.markers[i].points.size() > 0)
        freeNodesVis.markers[i].action = visualization_msgs::Marker::ADD;
      else
        freeNodesVis.markers[i].action = visualization_msgs::Marker::DELETE;
    }

    m_fmarkerPub.publish(freeNodesVis);
  }


  // finish pointcloud:
  if (publishPointCloud){
    sensor_msgs::PointCloud2 cloud;
    pcl::toROSMsg (pclCloud, cloud);
    cloud.header.frame_id = m_worldFrameId;
    cloud.header.stamp = rostime;
    m_pointCloudPub.publish(cloud);
  }

  if (publishBinaryMap)
    publishBinaryOctoMap(rostime);

  if (publishFullMap)
    publishFullOctoMap(rostime);


  double total_elapsed = (ros::WallTime::now() - startTime).toSec();
  ROS_DEBUG("Map publishing in TextureOctomapServer took %f sec", total_elapsed);

}


bool TextureOctomapServer::octomapBinarySrv(OctomapSrv::Request  &req,
                                    OctomapSrv::Response &res)
{
  ROS_INFO("Sending binary map data on service request");
  res.map.header.frame_id = m_worldFrameId;
  res.map.header.stamp = ros::Time::now();
  if (!octomap_msgs::binaryMapToMsg(*m_octree, res.map))
    return false;

  return true;
}

bool TextureOctomapServer::octomapFullSrv(OctomapSrv::Request  &req,
                                    OctomapSrv::Response &res)
{
  ROS_INFO("Sending full map data on service request");
  res.map.header.frame_id = m_worldFrameId;
  res.map.header.stamp = ros::Time::now();


  if (!octomap_msgs::fullMapToMsg(*m_octree, res.map))
    return false;

  return true;
}

bool TextureOctomapServer::clearBBXSrv(BBXSrv::Request& req, BBXSrv::Response& resp){
  point3d min = pointMsgToOctomap(req.min);
  point3d max = pointMsgToOctomap(req.max);

  for(OcTreeT::leaf_bbx_iterator it = m_octree->begin_leafs_bbx(min,max),
      end=m_octree->end_leafs_bbx(); it!= end; ++it){

    it->setLogOdds(octomap::logodds(m_thresMin));
    //			m_octree->updateNode(it.getKey(), -6.0f);
  }
  // TODO: eval which is faster (setLogOdds+updateInner or updateNode)
  m_octree->updateInnerOccupancy();

  publishAll(ros::Time::now());

  return true;
}

bool TextureOctomapServer::resetSrv(std_srvs::Empty::Request& req, std_srvs::Empty::Response& resp) {
  visualization_msgs::MarkerArray occupiedNodesVis;
  occupiedNodesVis.markers.resize(m_treeDepth +1);
  ros::Time rostime = ros::Time::now();
  m_octree->clear();
  // clear 2D map:
  m_gridmap.data.clear();
  m_gridmap.info.height = 0.0;
  m_gridmap.info.width = 0.0;
  m_gridmap.info.resolution = 0.0;
  m_gridmap.info.origin.position.x = 0.0;
  m_gridmap.info.origin.position.y = 0.0;

  ROS_INFO("Cleared octomap");
  publishAll(rostime);

  publishBinaryOctoMap(rostime);
  for (unsigned i= 0; i < occupiedNodesVis.markers.size(); ++i){

    occupiedNodesVis.markers[i].header.frame_id = m_worldFrameId;
    occupiedNodesVis.markers[i].header.stamp = rostime;
    occupiedNodesVis.markers[i].ns = "map";
    occupiedNodesVis.markers[i].id = i;
    occupiedNodesVis.markers[i].type = visualization_msgs::Marker::CUBE_LIST;
    occupiedNodesVis.markers[i].action = visualization_msgs::Marker::DELETE;
  }


  m_markerPub.publish(occupiedNodesVis);


  visualization_msgs::MarkerArray freeNodesVis;
  freeNodesVis.markers.resize(m_treeDepth +1);

  for (unsigned i= 0; i < freeNodesVis.markers.size(); ++i){

    freeNodesVis.markers[i].header.frame_id = m_worldFrameId;
    freeNodesVis.markers[i].header.stamp = rostime;
    freeNodesVis.markers[i].ns = "map";
    freeNodesVis.markers[i].id = i;
    freeNodesVis.markers[i].type = visualization_msgs::Marker::CUBE_LIST;
    freeNodesVis.markers[i].action = visualization_msgs::Marker::DELETE;
  }
  m_fmarkerPub.publish(freeNodesVis);

  return true;
}


void TextureOctomapServer::synthesizeView(const point3d& pos, const octomath::Quaternion& orient, 
                                          const unsigned int& h, const unsigned int& w,
                                          const float& fx, const float& fy,
                                          const float& cx, const float& cy,
                                          cv::Mat& image, cv::Mat& depth)
{
  image.create(h,w,CV_8U);
  depth.create(h,w,CV_32F);

  octomath::Quaternion imageToCam(0.5,-0.5,0.5,-0.5);
  octomath::Quaternion r = orient.inv() * imageToCam;

/*
#ifdef _OPENMP
  omp_set_num_threads(m_octree->keyrays.size());
  #pragma omp parallel for
#endif
*/

  // Cast Ray for each pixel
  for (unsigned v = 0; v < h; ++v)
  {
    for (unsigned u = 0; u < w; ++u)
    {
      float z = 1.0; // Use unit z for direction
      float x = (z/fx)*((u+0.5) - cx);
      float y = (z/fy)*((v+0.5) - cy);

      octomath::Vector3 pixDir(x,y,z);
      pixDir.normalize();
      octomath::Vector3 ray = r.rotate(pixDir);
      point3d end;
      
      // Cast ray to find voxel intersection
      bool ignoreUnknown = true;  
      bool found = m_octree->castRay(pos,ray,end,ignoreUnknown,m_maxRange);
      
      if(found)
      {
        // Set depth from distance to voxel center
        // Transform endpoint back to image coords
        octomath::Vector3 p = end - pos;
        p = r.inv().rotate(p);
        depth.at<float>(v,u) = p.z();

        // Compute intensity of voxel face that ray intersects
        image.at<unsigned char>(v,u) = m_octree->getTexturePoint(end,ray,pos); 
      }
      else
      {
        image.at<unsigned char>(v,u) = 0;
        depth.at<float>(v,u) = 0.0;
      }
    }
  }

  return;
}

bool TextureOctomapServer::synthesizeViewsSrv(TextureSrv::Request &req, TextureSrv::Response &res) {

  ros::Time start = ros::Time::now();
  std::vector<sensor_msgs::Image> images;
  std::vector<sensor_msgs::Image> depths;

  for(unsigned i = 0; i < req.poses.size(); i++)
  {
    // Convert point
    point3d pos = pointMsgToOctomap(req.poses.at(i).position);
    // Convert orientation
    octomath::Quaternion orient(req.poses.at(i).orientation.w, 
                                req.poses.at(i).orientation.x, 
                                req.poses.at(i).orientation.y,
                                req.poses.at(i).orientation.z);

    // Create some opencv mats to hold the output
    cv::Mat im, d;

    // Call the view synthesizer function for this pose
    synthesizeView(pos,orient,req.h,req.w,req.fx,req.fy,req.cx,req.cy,im,d);

    // Pack up the output 
    cv_bridge::CvImagePtr im_ptr(new cv_bridge::CvImage);
    im_ptr->image = im;
    im_ptr->encoding = "mono8";
    cv_bridge::CvImagePtr d_ptr(new cv_bridge::CvImage);
    d_ptr->image = d;
    d_ptr->encoding = sensor_msgs::image_encodings::TYPE_32FC1;

    sensor_msgs::Image image = *(im_ptr->toImageMsg());
    images.push_back(image);
    sensor_msgs::Image depth = *(d_ptr->toImageMsg());
    depths.push_back(depth);
  }
  res.images = images;
  res.depths = depths;
  
  ros::Time end = ros::Time::now();
  double elapsed = (end-start).toSec();
  ROS_INFO("Time to synthesize %d views: %f", (int) req.poses.size(), elapsed);
  
  return true;
}

void TextureOctomapServer::publishBinaryOctoMap(const ros::Time& rostime) const{

  Octomap map;
  map.header.frame_id = m_worldFrameId;
  map.header.stamp = rostime;

  if (octomap_msgs::binaryMapToMsg(*m_octree, map))
    m_binaryMapPub.publish(map);
  else
    ROS_ERROR("Error serializing OctoMap");
}

void TextureOctomapServer::publishFullOctoMap(const ros::Time& rostime) const{

  Octomap map;
  map.header.frame_id = m_worldFrameId;
  map.header.stamp = rostime;

  if (octomap_msgs::fullMapToMsg(*m_octree, map))
    m_fullMapPub.publish(map);
  else
    ROS_ERROR("Error serializing OctoMap");

}

void TextureOctomapServer::handlePreNodeTraversal(const ros::Time& rostime){
  if (m_publish2DMap){
    // init projected 2D map:
    m_gridmap.header.frame_id = m_worldFrameId;
    m_gridmap.header.stamp = rostime;
    nav_msgs::MapMetaData oldMapInfo = m_gridmap.info;

    // TODO: move most of this stuff into c'tor and init map only once (adjust if size changes)
    double minX, minY, minZ, maxX, maxY, maxZ;
    m_octree->getMetricMin(minX, minY, minZ);
    m_octree->getMetricMax(maxX, maxY, maxZ);

    octomap::point3d minPt(minX, minY, minZ);
    octomap::point3d maxPt(maxX, maxY, maxZ);
    octomap::OcTreeKey minKey = m_octree->coordToKey(minPt, m_maxTreeDepth);
    octomap::OcTreeKey maxKey = m_octree->coordToKey(maxPt, m_maxTreeDepth);
    
    ROS_DEBUG("MinKey: %d %d %d / MaxKey: %d %d %d", minKey[0], minKey[1], minKey[2], maxKey[0], maxKey[1], maxKey[2]);

    // add padding if requested (= new min/maxPts in x&y):
    double halfPaddedX = 0.5*m_minSizeX;
    double halfPaddedY = 0.5*m_minSizeY;
    minX = std::min(minX, -halfPaddedX);
    maxX = std::max(maxX, halfPaddedX);
    minY = std::min(minY, -halfPaddedY);
    maxY = std::max(maxY, halfPaddedY);
    minPt = octomap::point3d(minX, minY, minZ);
    maxPt = octomap::point3d(maxX, maxY, maxZ);

    OcTreeKey paddedMaxKey;
    if (!m_octree->coordToKeyChecked(minPt, m_maxTreeDepth, m_paddedMinKey)){
      ROS_ERROR("Could not create padded min OcTree key at %f %f %f", minPt.x(), minPt.y(), minPt.z());
      return;
    }
    if (!m_octree->coordToKeyChecked(maxPt, m_maxTreeDepth, paddedMaxKey)){
      ROS_ERROR("Could not create padded max OcTree key at %f %f %f", maxPt.x(), maxPt.y(), maxPt.z());
      return;
    }

    ROS_DEBUG("Padded MinKey: %d %d %d / padded MaxKey: %d %d %d", m_paddedMinKey[0], m_paddedMinKey[1], m_paddedMinKey[2], paddedMaxKey[0], paddedMaxKey[1], paddedMaxKey[2]);
    assert(paddedMaxKey[0] >= maxKey[0] && paddedMaxKey[1] >= maxKey[1]);

    m_multires2DScale = 1 << (m_treeDepth - m_maxTreeDepth);
    m_gridmap.info.width = (paddedMaxKey[0] - m_paddedMinKey[0])/m_multires2DScale +1;
    m_gridmap.info.height = (paddedMaxKey[1] - m_paddedMinKey[1])/m_multires2DScale +1;

    int mapOriginX = minKey[0] - m_paddedMinKey[0];
    int mapOriginY = minKey[1] - m_paddedMinKey[1];
    assert(mapOriginX >= 0 && mapOriginY >= 0);

    // might not exactly be min / max of octree:
    octomap::point3d origin = m_octree->keyToCoord(m_paddedMinKey, m_treeDepth);
    double gridRes = m_octree->getNodeSize(m_maxTreeDepth);
    m_projectCompleteMap = (!m_incrementalUpdate || (std::abs(gridRes-m_gridmap.info.resolution) > 1e-6));
    m_gridmap.info.resolution = gridRes;
    m_gridmap.info.origin.position.x = origin.x() - gridRes*0.5;
    m_gridmap.info.origin.position.y = origin.y() - gridRes*0.5;
    if (m_maxTreeDepth != m_treeDepth){
      m_gridmap.info.origin.position.x -= m_res/2.0;
      m_gridmap.info.origin.position.y -= m_res/2.0;
    }
    
    // workaround for  multires. projection not working properly for inner nodes:
    // force re-building complete map
    if (m_maxTreeDepth < m_treeDepth)
      m_projectCompleteMap = true;


    if(m_projectCompleteMap){
      ROS_DEBUG("Rebuilding complete 2D map");
      m_gridmap.data.clear();
      // init to unknown:
      m_gridmap.data.resize(m_gridmap.info.width * m_gridmap.info.height, -1);

    } else {

       if (mapChanged(oldMapInfo, m_gridmap.info)){
          ROS_DEBUG("2D grid map size changed to %dx%d", m_gridmap.info.width, m_gridmap.info.height);
          adjustMapData(m_gridmap, oldMapInfo);
       }
       nav_msgs::OccupancyGrid::_data_type::iterator startIt;
       size_t mapUpdateBBXMinX = std::max(0, (int(m_updateBBXMin[0]) - int(m_paddedMinKey[0]))/int(m_multires2DScale));
       size_t mapUpdateBBXMinY = std::max(0, (int(m_updateBBXMin[1]) - int(m_paddedMinKey[1]))/int(m_multires2DScale));
       size_t mapUpdateBBXMaxX = std::min(int(m_gridmap.info.width-1), (int(m_updateBBXMax[0]) - int(m_paddedMinKey[0]))/int(m_multires2DScale));
       size_t mapUpdateBBXMaxY = std::min(int(m_gridmap.info.height-1), (int(m_updateBBXMax[1]) - int(m_paddedMinKey[1]))/int(m_multires2DScale));

       assert(mapUpdateBBXMaxX > mapUpdateBBXMinX);
       assert(mapUpdateBBXMaxY > mapUpdateBBXMinY);

       size_t numCols = mapUpdateBBXMaxX-mapUpdateBBXMinX +1;

       // test for max idx:
       uint max_idx = m_gridmap.info.width*mapUpdateBBXMaxY + mapUpdateBBXMaxX;
       if (max_idx  >= m_gridmap.data.size())
         ROS_ERROR("BBX index not valid: %d (max index %zu for size %d x %d) update-BBX is: [%zu %zu]-[%zu %zu]", max_idx, m_gridmap.data.size(), m_gridmap.info.width, m_gridmap.info.height, mapUpdateBBXMinX, mapUpdateBBXMinY, mapUpdateBBXMaxX, mapUpdateBBXMaxY);

       // reset proj. 2D map in bounding box:
       for (unsigned int j = mapUpdateBBXMinY; j <= mapUpdateBBXMaxY; ++j){
          std::fill_n(m_gridmap.data.begin() + m_gridmap.info.width*j+mapUpdateBBXMinX,
                      numCols, -1);
       }

    }

  }

}

void TextureOctomapServer::handlePostNodeTraversal(const ros::Time& rostime){

  if (m_publish2DMap)
    m_mapPub.publish(m_gridmap);
}

void TextureOctomapServer::handleOccupiedNode(const OcTreeT::iterator& it){

  if (m_publish2DMap && m_projectCompleteMap){
    update2DMap(it, true);
  }
}

void TextureOctomapServer::handleFreeNode(const OcTreeT::iterator& it){

  if (m_publish2DMap && m_projectCompleteMap){
    update2DMap(it, false);
  }
}

void TextureOctomapServer::handleOccupiedNodeInBBX(const OcTreeT::iterator& it){

  if (m_publish2DMap && !m_projectCompleteMap){
    update2DMap(it, true);
  }
}

void TextureOctomapServer::handleFreeNodeInBBX(const OcTreeT::iterator& it){

  if (m_publish2DMap && !m_projectCompleteMap){
    update2DMap(it, false);
  }
}

void TextureOctomapServer::update2DMap(const OcTreeT::iterator& it, bool occupied){

  // update 2D map (occupied always overrides):

  if (it.getDepth() == m_maxTreeDepth){
    unsigned idx = mapIdx(it.getKey());
    if (occupied)
      m_gridmap.data[mapIdx(it.getKey())] = 100;
    else if (m_gridmap.data[idx] == -1){
      m_gridmap.data[idx] = 0;
    }

  } else{
    int intSize = 1 << (m_maxTreeDepth - it.getDepth());
    octomap::OcTreeKey minKey=it.getIndexKey();
    for(int dx=0; dx < intSize; dx++){
      int i = (minKey[0]+dx - m_paddedMinKey[0])/m_multires2DScale;
      for(int dy=0; dy < intSize; dy++){
        unsigned idx = mapIdx(i, (minKey[1]+dy - m_paddedMinKey[1])/m_multires2DScale);
        if (occupied)
          m_gridmap.data[idx] = 100;
        else if (m_gridmap.data[idx] == -1){
          m_gridmap.data[idx] = 0;
        }
      }
    }
  }


}



bool TextureOctomapServer::isSpeckleNode(const OcTreeKey&nKey) const {
  OcTreeKey key;
  bool neighborFound = false;
  for (key[2] = nKey[2] - 1; !neighborFound && key[2] <= nKey[2] + 1; ++key[2]){
    for (key[1] = nKey[1] - 1; !neighborFound && key[1] <= nKey[1] + 1; ++key[1]){
      for (key[0] = nKey[0] - 1; !neighborFound && key[0] <= nKey[0] + 1; ++key[0]){
        if (key != nKey){
          OcTreeNode* node = m_octree->search(key);
          if (node && m_octree->isNodeOccupied(node)){
            // we have a neighbor => break!
            neighborFound = true;
          }
        }
      }
    }
  }

  return neighborFound;
}

void TextureOctomapServer::reconfigureCallback(octomap_server::OctomapServerConfig& config, uint32_t level){
  if (m_maxTreeDepth != unsigned(config.max_depth)){
    m_maxTreeDepth = unsigned(config.max_depth);

    publishAll();
  }


}

void TextureOctomapServer::adjustMapData(nav_msgs::OccupancyGrid& map, const nav_msgs::MapMetaData& oldMapInfo) const{
  if (map.info.resolution != oldMapInfo.resolution){
    ROS_ERROR("Resolution of map changed, cannot be adjusted");
    return;
  }

  int i_off = int((oldMapInfo.origin.position.x - map.info.origin.position.x)/map.info.resolution +0.5);
  int j_off = int((oldMapInfo.origin.position.y - map.info.origin.position.y)/map.info.resolution +0.5);

  if (i_off < 0 || j_off < 0
      || oldMapInfo.width  + i_off > map.info.width
      || oldMapInfo.height + j_off > map.info.height)
  {
    ROS_ERROR("New 2D map does not contain old map area, this case is not implemented");
    return;
  }

  nav_msgs::OccupancyGrid::_data_type oldMapData = map.data;

  map.data.clear();
  // init to unknown:
  map.data.resize(map.info.width * map.info.height, -1);

  nav_msgs::OccupancyGrid::_data_type::iterator fromStart, fromEnd, toStart;

  for (int j =0; j < int(oldMapInfo.height); ++j ){
    // copy chunks, row by row:
    fromStart = oldMapData.begin() + j*oldMapInfo.width;
    fromEnd = fromStart + oldMapInfo.width;
    toStart = map.data.begin() + ((j+j_off)*m_gridmap.info.width + i_off);
    copy(fromStart, fromEnd, toStart);

//    for (int i =0; i < int(oldMapInfo.width); ++i){
//      map.data[m_gridmap.info.width*(j+j_off) +i+i_off] = oldMapData[oldMapInfo.width*j +i];
//    }

  }

}


std_msgs::ColorRGBA TextureOctomapServer::heightMapColor(double h) {

  std_msgs::ColorRGBA color;
  color.a = 1.0;
  // blend over HSV-values (more colors)

  double s = 1.0;
  double v = 1.0;

  h -= floor(h);
  h *= 6;
  int i;
  double m, n, f;

  i = floor(h);
  f = h - i;
  if (!(i & 1))
    f = 1 - f; // if i is even
  m = v * (1 - s);
  n = v * (1 - s * f);

  switch (i) {
    case 6:
    case 0:
      color.r = v; color.g = n; color.b = m;
      break;
    case 1:
      color.r = n; color.g = v; color.b = m;
      break;
    case 2:
      color.r = m; color.g = v; color.b = n;
      break;
    case 3:
      color.r = m; color.g = n; color.b = v;
      break;
    case 4:
      color.r = n; color.g = m; color.b = v;
      break;
    case 5:
      color.r = v; color.g = m; color.b = n;
      break;
    default:
      color.r = 1; color.g = 0.5; color.b = 0.5;
      break;
  }

  return color;
}


void TextureOctomapServer::trackChanges() {
  KeyBoolMap::const_iterator startPnt = m_octree->changedKeysBegin();
  KeyBoolMap::const_iterator endPnt = m_octree->changedKeysEnd();

  pcl::PointCloud<pcl::PointXYZI> changedCells = pcl::PointCloud<pcl::PointXYZI>();

  int c = 0;
  for (KeyBoolMap::const_iterator iter = startPnt; iter != endPnt; ++iter) {
    c++;
    OcTreeNode* node = m_octree->search(iter->first);

    bool occupied = m_octree->isNodeOccupied(node);

    pcl::PointXYZI pnt;
    pnt.x = m_octree->keyToCoord(iter->first.k[0]);
    pnt.y = m_octree->keyToCoord(iter->first.k[1]);
    pnt.z = m_octree->keyToCoord(iter->first.k[2]);

    if (occupied) {
      pnt.intensity = 1000;
    }
    else {
      pnt.intensity = -1000;
    }

    changedCells.push_back(pnt);
  }

  sensor_msgs::PointCloud2 changed;
  pcl::toROSMsg(changedCells, changed);
  changed.header.frame_id = "/talker/changes";
  changed.header.stamp = ros::Time().now();
  pubChangeSet.publish(changed);
  ROS_DEBUG("[server] sending %d changed entries", (int)changedCells.size());

  m_octree->resetChangeDetection();
  ROS_DEBUG("[server] octomap size after updating: %d", (int)m_octree->calcNumNodes());
}





}



