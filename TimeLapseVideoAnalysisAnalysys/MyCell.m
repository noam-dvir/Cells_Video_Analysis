classdef MyCell < handle
   properties
      lbl
      frame
      area
      perimeter
      centroid(1,2)
      prevCentroid(1,2)
      driftDistance
      driftAngle
      isJustCreated
      neighborsList(2,10)
      % borderList
   end
   methods
      function obj = MyCell()

      end
      
      function m = toMatrix(obj)
            m = [obj.frame obj.lbl obj.area obj.perimeter obj.centroid obj.prevCentroid obj.driftDistance obj.driftAngle obj.isJustCreated];
            m = [m obj.neighborsList(1,:) obj.neighborsList(2,:)];
      end
      
      function m = fromMatrix(obj, line)
            obj.frame = line(1);
            obj.lbl = line(2);
            obj.area = line(3);
            obj.perimeter = line(4);
            obj.centroid = [line(5:6)];
            obj.prevCentroid = [line(7:8)];
            obj.driftDistance = line(9);
            obj.driftAngle= line(10);
            obj.isJustCreated= line(11);
            obj.neighborsList = [line(12:21) ; line(22:end)];
      end
      
      function y = getDriftCord(obj)
          x1 = obj.centroid(1);
          y1 = obj.centroid(2);
          x2 = obj.prevCentroid(1);
          y2 = obj.prevCentroid(2);
          
          y = [x1-x2, y1-y2]; 
      end
      
      function y = getDriftDistance(obj)
          if (obj.driftAngle==-1)
              delta = getDriftCord(obj);
              obj.driftDistance = sqrt(delta(1).^2 + delta(2).^2);
          end
          y = obj.driftDistance;
      end
      
      function y = getDriftAngle(obj)
          if (obj.driftAngle==-1)
              delta = getDriftCord(obj);
              obj.driftAngle = atan2d(delta(2),delta(1));
          end
          y = obj.driftAngle;
      end
      
      function y = getNeighborsSum(obj)
          y = sum(obj.neighborsList(1,:) > 0);
      end
      
      function y = getNeighborsList(obj)
          y = obj.neighborsList(1,:);
      end
      
      function y = getNeighborsEdgeSize(obj, lbl)
          y = 0;
          neighborIndex = find(getNeighborsList(obj)==lbl)
          if neighborIndex
              y = obj.neighborsList(2, neighborIndex);
          end
      end
      
   end
end